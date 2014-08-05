#################################################################################
#
#       PROJECT VISUFRONT - Read validated ISIIS data and plot them
#
#            Created by Robin Faillettaz on date 13/07/2014
#     UPMC - Laboratoire d'Océanograhie de Villefranche-sur-Mer (LOV)       
#
#################################################################################

# Notes
 
# A .Rprofile file should be set personnaly to the directory where the data are stored (dropbox)
# for CC5, casts are inverted (01 is offshore and 54 is ashore) and have to be reverted with nbCast+1 - castNb




# For now locate data from the dropbox repo manually
dir <- "/Users/faillettaz/Dropbox/visufront-data/"


# Load libraries
library("plyr")
library("ggplot2")
library("stringr")
library("reshape2")
library("grid")
library("gridExtra")
library("akima")
library("fields")
library("oce")   


source("lib_zooprocess.R")
source("lib_plot_rf.R")
source("lib_process.R")


# set options for ploting
opts <- theme(axis.title.y= element_text(angle=90, vjust=0.5, size=18), axis.title.x= element_text(angle=0, vjust=0.5, size=18)) + 
theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=15), axis.text.y  = element_text(angle=0, vjust=0.5, size=15))



# --------------------------------------------------------------------
#       READ ZOOPROCESS FILES
# --------------------------------------------------------------------
 {


# list files
files <- list.files(str_c(dir, "zooprocess/"), full = T)


# Select dat1 files
# dat1 <- files[which(str_detect(files, "_dat1.txt") == T)] # For all files
dat1 <- files[which(str_detect(files, "byR_dat1.txt") == T)] # for certain files only

# select datfiles
datfiles <- files[which(str_detect(files, "_datfile.txt") == T)]


# read selected dat1 files 
pids <- adply(dat1, 1, function(x) {

pid <- read.pid(x)

# make validation column name uniform
names <- names(pid)

#names <- str_replace_all(names, "Valid", "Valid")
# force latest validation to have the name "Valid"
validationColumns <- which(str_detect(names, "Valid"))
if (length(validationColumns) == 1) {
    names[validationColumns] <- "Valid"
}
if (length(validationColumns) > 1) {
    names[max(validationColumns)] <- "Valid"    
}

names(pid) <- names
return(pid)

}, .progress="text")  
pids <- pids[, -1]


 }


# --------------------------------------------------------------------
#       CLEAN NAMES AND MERGE GROUPS
# --------------------------------------------------------------------
 {


# Check "not_found" category
pids[which(pids$Valid == 'not_found'), ]

# correspond to images lost during sorting, get rid of them
if (length(pids[which(pids$Valid == 'not_found'), "Valid"] > 0)) {
    message("Removing lost images (not found)")
    pids <- pids[-which(pids$Valid == 'not_found'), ]
}



# Merge groups into large taxonomical groups
# --------------------------------------------------------------------

    
    
pids$groups <- pids$Valid 

# Merge appendicularians
pids[which(str_detect(pids$Valid, "append") == T), "groups"] <- "appendicularians"

# merge radiolarians names
pids[which(str_detect(pids$Valid, "radiolarian") == T), "groups"] <- "radiolarians"

# merge fish_like
pids[which(str_detect(pids$Valid, "fish") == T), "groups"] <- "fish_like"

# Merde crustaceans
pids[which(str_detect(pids$Valid, "crust") == T), "groups"] <- "crustaceans"

# Merde jellyfishes
pids[which(str_detect(pids$Valid, "jelly") == T | pids$Valid %in% "ephyrae"), "groups"] <- "jellyfishes"

# Merde pteropods
pids[which(str_detect(pids$Valid, "pterop") == T), "groups"] <- "pteropods"

# Merde aggregates
pids[which(str_detect(pids$Valid, "det") == T), "groups"] <- "detritus"

# Merde siphonophores
pids[which(str_detect(pids$Valid, "siphos") == T), "groups"] <- "siphonophores"

# Merde others
pids[which(pids$Valid %in% c("bad_focus", "duplicates", "unidentified", "unidentified_of_interest")), "groups"] <- "others"

sort(unique(pids$groups))


 }


#---------------------------------------------------------------------
#       EXTRACT BIOLOGICAL INFO
#---------------------------------------------------------------------
 {


# set depth BIN 
binsize <- 1
pids$DepthBin <- round_any(pids$Depth, binsize)

# get raw abundance
bio <- ddply(pids, ~Valid+Label+DepthBin, function(x) {
    sum(na.omit(x$Valid==paste(x$Valid)))
    })
bio <- rename(bio, c("V1" = "Abund"))
head(bio)

# check for NAs
if (dim(bio[which(is.na(bio)), ])[1] > 0) {
    message("Removing NAs")
    bio <- bio[-which(is.na(bio)), ]
    }

head(bio)



# join with depth bin to have a line for each depth
# Set all possibilities of depth, valid and label and add corresponding abundances
grid <- expand.grid(sort(unique(bio$Valid)), sort(unique(bio$Label)), sort(unique(bio$DepthBin)))
head(grid)
names(grid) <- c("Valid", "Label", "DepthBin")

# join the grid with abundances
bioFull <- join(grid, bio)


# replace NAs by 0
bioFull[which(is.na(bioFull$Abund)), "Abund"] <- 0 
length(which(is.na(bioFull)))  # if 0 --> OK


# get number of a certain group
sum(bioFull[which(bioFull$Valid == "det_aggregates"), "Abund"])


# get cast number from the profil name (vignette label)
bioFull$cast <- 54 - (as.numeric(str_split_fixed(bioFull$Label, fixed("_"), 3)[, 3]) * 2 - 1)


# plot everything
ggplot(bioFull[which(bioFull$Valid %in% c("fish_larvae", "fish_like", "chaetognaths", "siphos_calycophore")), ]) + 
geom_point(aes(x = Abund, y = -DepthBin, colour = Label, group = Label)) + 
#geom_path(aes(x = Abund, y = -DepthBin, colour = Label)) + 
facet_grid(.~Valid, scales="free_x")


 }


#---------------------------------------------------------------------
#       COMPUTE VOLUME SAMPLED PER DEPTH BIN
#---------------------------------------------------------------------
 {
     
vol <- adply(datfiles[14:15], 1, function(x) {
    # read datfiles
    t <- read.table(x, sep="\t", h=F)

    # get depth
    t$depth <- as.numeric(substring(t$V3, 1, 5))/10
    class(t$depth)

    # round depth to bin per m
    t$DepthBin <- round_any(t$depth, binsize)

    # get number of images per depth bin
    img <- ddply(t, ~DepthBin, function(x){
        l <- length(x$V1)
    })
    img <- rename(img, c('V1' = 'nb.img'))

    
    # keep only uniques of each line
    #vol <- v[!duplicated(v[, c("DepthBin", "nb")]), c("DepthBin", "nb")]

    img$vol.m3 <- img$nb * 7.78 / 1000
    img$cast <- 54 - (as.numeric(str_split_fixed(x, fixed("_"), 4)[, 3]) * 2 -1) 
    # 55 - should be removed when processing other transects
    
    
    return(img)
    
}, .progress="text")
    vol <- vol[, -1]
    
head(vol)
length(which(is.na(vol$vol.m3)))  # if 0 --> OK

 }


#---------------------------------------------------------------------
#       COMPUTE ABUNDANCES (VOLUMES)
#---------------------------------------------------------------------
 {
 
head(bioFull)
head(vol)
bioFull <- join(bioFull, vol)
bioFull$abund.m3 <- bioFull$Abund / bioFull$vol.m3


# NAs
length(which(is.na(bioFull$nb.img)))  # if 0 --> OK

if (dim(bioFull[which(is.na(bioFull$abund.m3)), ])[1] > 0) {
    message("removing NAs")
    bioFull <- bioFull[-which(is.na(bioFull$nb.img)), ]   # Only some column work, weird. It seems to be from the vol object after joining
}

 }
 
 
#---------------------------------------------------------------------
#       READ PHYSICAL DATA
#---------------------------------------------------------------------
 {
 

# Read physical data from transect 4
phy <- read.csv(str_c(dir, "ISIIShydro/transects/cross_current_5/isiis.csv"), header=T, sep=",")
head(phy)
# delete first line (data from previous transect)
phy <- phy[-1, ]

# Change names to same number of character for each (for plots)
phy <- rename(phy, c("Temp.C" = "Temp.celsius"))

# select only data from upcasts --> Keep first down cast below 25m to improve interpolation
data <- phy[which(phy$down.up %in% "up" | phy$down.up %in% "down" & phy$Depth.m > 28), ]


# check if seems ok
ggplot(data, aes(x=distanceFromVlfr, y=-Depth.m, colour=Salinity.PPT))  + geom_point() # Yes !

 }


#---------------------------------------------------------------------
#        GLIDER DATA ---- BROKEN ------
#---------------------------------------------------------------------
 {
# # Add glider data that to complete the profile
# g <- read.csv("~/Desktop/PhD/VISUFRONT/ISIIS/glider/glier-d25.csv", sep=",", header=T)
# head(g)
# 
# 
# # check path
# ggplot(g) + geom_path(aes(x=distance, y=ipressure, colour= substring(dateTime, 1, 10)))
# 
# # get data from cast 1 to 4
# g[which(round(g$distance) == 7), ]
# min(g$distance)
# 
# 
# # select data to keep 
# # keep only variable of interest
# g <- g[which(g$ipressure < 103 & g$distance < 7.2), names(g) %in% c("distance", "ipressure", "salinity")]
# 
# 
# d_short <- d[ , names(d) %in% c("distanceFromVlfr", "Depth.m", "Salinity.PPT")]
# head(d_short)
# 
# colnames(g) <- names(d_short)
# head(g)
# 
# dg <- rbind(d_short, g)
# head(dg)
# 
# 
# # compute interpolation 
# # interpolate all variables
# dm <- melt(dg, id.vars=c("Depth.m", "distanceFromVlfr"), measure.vars=c("Salinity.PPT"))#, "Temp.C", "Fluoro.volts", "Oxygen.ml.l"))
# 
# 
# # First interpolation over a large grid
# di <- ddply(dm, ~variable, function(x) {
#     x <- na.omit(x)
#     xi <- interp.dist(x=x$distanceFromVlfr, y=x$Depth.m, z=x$value, duplicate="mean", x.step=500, y.step=2.5, anisotropy=1300)
# }, .progress="text")
# 
# di <- rename(di, c("x"="distance", "y"="Depth.m"))
# 
# plots <- dlply(di, ~variable, function(x) {
#         ggplot(x, aes(x=distance, y=-Depth.m)) +
#         #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
#         geom_tile(aes(fill=value), na.rm=T) +
#         stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
#         scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
#         scale_x_continuous(expand=c(0,0)) +
#         scale_y_continuous(expand=c(0,0))
#         })
# 
# do.call(grid.arrange, c(plots,list(ncol=1)))

 }

#---------------------------------------------------------------------
#        PROCESS PHYSICAL DATA IF NO GLIDER DATA 
#---------------------------------------------------------------------
 {


# compute interpolation 
# interpolate all variables
dm <- melt(data, id.vars=c("Depth.m", "distanceFromVlfr"), measure.vars=c("Salinity.PPT", "Temp.celsius", "Fluoro.volts", "Oxygen.ml.l"))



# First interpolation over a large grid
di <- ddply(dm, ~variable, function(x) {
    x <- na.omit(x)
    xi <- interp.dist(x=x$distanceFromVlfr, y=x$Depth.m, z=x$value, duplicate="mean", x.step=500, y.step=2.5, anisotropy=1300)
}, .progress="text")

di <- rename(di, c("x"="distance", "y"="Depth.m"))

plots <- dlply(di, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=-Depth.m)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_tile(aes(fill=value), na.rm=T) +
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
        })

do.call(grid.arrange, c(plots,list(ncol=1)))




# Second interpolation over a finer grid for ploting
#----------------------------------------------------

# if updates in the lib_plot
#source("data/lib_plot.R")


# Second interpolation for all variables
di2 <- ddply(di, ~variable, function(x) {
    x <- na.omit(x)
    xi <- interp.smooth(x=x$distance, y=x$Depth.m, z=x$value, x.step = 0.1, y.step = 0.1)
}, .progress="text")


di2 <- rename(di2, c("x"="distance", "y"="Depth.m"))

plots <- dlply(di2, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=Depth.m)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(x$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("Distance from shore", expand=c(0,0)) +
        scale_y_reverse("Depth (m)", expand=c(0,0)) +
        opts
        })

do.call(grid.arrange, c(plots,list(ncol=1)))


 }


# --------------------------------------------------------------------
#       MERGE PHYSICAL AND BIOLOGICAL DATA
# --------------------------------------------------------------------
 
# Merge interpolated physical and biological data
# -------------------------------------------------------------------- 
 {

# NB : Bio data need to be localised in x, y, z --> I will bin the depth by profile
# and get the corresponding x and y

head(bioFull)
any(is.na(bioFull))


# NB : Bio data need to be localised in x, y, z --> I will bin the depth by profile
# and get the corresponding x and y




# !!!!!!!!!!!!!!!!!!  NEED TO BE FIXED !!!!!!!!!!!!!!!!!!!!!!
# cast are inverted for dat1 files and physical files 









# select only physical data from the profil 05 (cast 9)
biophy <- ddply(bioFull, ~cast, function(x) {
    
    # x <- bioFull[which(bioFull$cast == 1), ]
    
    phyx <- phy[which(phy$cast == unique(x$cast)), ]
    
    # Get biological data
    biox <- bioFull[which(bioFull$cast == unique(x$cast)), ]

    # cast the df
    dC <- dcast(biox, DepthBin ~ Valid, value.var = "abund.m3")
    length(which(is.na(dC)))
    
    
    
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # ADD A DCAST FOR RAW ABUND AND NOT ONLY ABUND.M3
    
    

    
    
    
    print(unique(phyx$cast))
   
    # Bin depth = select middle of each meter
    # set depth BIN 
    bin <- seq(from=1, to = ceiling(max(phyx$Depth.m)), by=binsize)
    
    # get index of the min of each point
    min <- adply(bin, 1, function(m){
        which.min(abs(phyx$Depth.m-m))
        })
        
    which(duplicated(min$V1))
    
    # select only the physical data corresponding to the center of each bin
    phyxMin <- phyx[min$V1, ]
    
    
    # IF DUPLICATEDS, it generates NAs at the end that are hard to get rid of
    if (length(which(duplicated(phyxMin))) > 0 ) {
        
        # for factor variables
        phyxMin[which(duplicated(phyxMin)), names(phyxMin) %in% c("dateTime", "down.up", "dateTimeMsec")] <- phyxMin[which(duplicated(phyxMin))-1, names(phyxMin) %in% c("dateTime", "down.up", "dateTimeMsec")]
        
        # for numerical variable, do the mean of the previous and next
        phyxMin[which(duplicated(phyxMin)), !names(phyxMin) %in% c("dateTime", "down.up", "dateTimeMsec")] <- 
        (phyxMin[which(duplicated(phyxMin))-1, !names(phyxMin) %in% c("dateTime", "down.up", "dateTimeMsec")] +  phyxMin[which(duplicated(phyxMin))+1, !names(phyxMin) %in% c("dateTime", "down.up", "dateTimeMsec")]) / 2
        
     }
     
     # net very good but force bin
     phyxMin$DepthBin <- bin 
     
     
#     # Shoud be done this way
#     phyxMin$DepthBin <-  round_any(phyxMin$Depth.m, binsize)
#
#     if (phyxMin$DepthBin[length(phyxMin$DepthBin)] == 102) {
#        phyxMin$DepthBin[length(phyxMin$DepthBin)] <- 103
#    }
    
    #dim(phyxMin)
    #dim(dC)
    
    #head(phyxMin)
    #head(dC)
    
    
    
    # join the 2
    biophyx <- join(dC, phyxMin)


    #biophyx[which(is.na(biophyx)), ]
    return(biophyx)
    
}, .progress="text")

# NB : If warnings it's because x$cast needs to be unique (solved)
# Some problem with some bin that doesn't exist in the physical data
# These data should be interpolated later on

head(biophy)
dim(biophy)
unique(biophy$cast)
table(biophy$DepthBin)
table(biophy$cast)

dim(biophy[which(is.na(biophy$Salinity.PPT)), ])


# remove NAs
if (dim(biophy[which(is.na(biophy$Salinity.PPT)), ])[1] > 0) {
    message("Removing NAs")
    biophy <- biophy[-which(is.na(biophy$Salinity.PPT)), ]   # other columns are <NA> so not recognized as NAs
}

head(biophy)


# try w/ 1 variable only
sal <- di2[which(di2$variable == "Salinity.PPT"), ]


# Melt biophy plot plots
biophyM <- melt(biophy, id.vars = c("DepthBin", "cast", "dateTimeMsec", "Pressure.dbar", "Depth.m", "Temp.celsius", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "Salinity.PPT", "Density", "dateTime", "lat", "lon", "distanceFromStart", "distanceFromVlfr",   "distanceFromShore", "down.up"), variable.name = "taxa", value.name = "abund.m3")
head(biophyM)
tail(biophyM)

# change NA value of abundance, due to absence of data on the grid, to 0 since it's absences 
biophyM[which(is.na(biophyM$abund.m3)), ]
dim(which(is.na(biophyM$abund.m3)))
biophyM$abund.m3[which(is.na(biophyM$abund.m3))] <- 0


# plot fish larvae
data <- biophyM[which(biophyM$taxa == "fish_larvae"), ]
ggplot() + 
geom_raster(aes(x=distance, y=-Depth.m, fill=value), data= sal , na.rm=T, ) +
stat_contour(aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE, data=sal) +
#geom_path(aes(x=distanceFromVlfr, y=-Depth.m, group=cast), size=0.6, data=phy) +
geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3, group=cast))+
scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
scale_x_continuous(expand=c(0,0)) +
scale_y_continuous(expand=c(0,0)) +
scale_size_area()

# IT WORKS, GREAT, but physical data should be extrapolated on the sides

 }

# Smooth and extrapolate physical data
# --------------------------------------------------------------------
 {
     
# Here, we use image.smooth that does the extrapolation on the left of the plot and make it nicer

head(biophy)

# interpolate like above for each variable
i <- dlply(di, ~variable, function(x) {
    # get rid of NAs
    x <- na.omit(x)
    
    # select x y and z
    x1 <- x$distance
    y1 <- x$Depth.m
    z1 <- x$value
    
    # run linear interpolation
    i <- interp(x=x1, y=y1, z=z1, xo=seq(0, max(x1), by=0.1), yo=seq(0, max(y1), by=0.1),  linear = T, duplicate = "mean", extrap=F)
    
    return(i)
}, .progress="text")


# Create the grid for the nex interpolation w/ extrapolation
grid <- data.frame(x=seq(0, max(i$Salinity.PPT$x), by=0.1), y=seq(0, max(i$Salinity.PPT$y), length=length(seq(0, max(i$Salinity.PPT$x), by=0.1))))

iS <- ldply(i, function(x) {
    # run smoothing
    smooth <- image.smooth(x, grid=grid, theta=0.28)
    
    # pass the list to df
    out <- melt(smooth$z, varnames=c("x","y"))
    out$x <- smooth$x[out$x]
    out$y <- smooth$y[out$y]
    out <- rename(out, c("x"="distance", "y"="Depth.m"))
    
    return(out)
    
}, .progress="text")

head(iS)
unique(iS$variable)

# save data for next time
#save(is, file = "physics.Rdata")

 }


# plot it 
# --------------------------------------------------------------------

 {

# short way but axis labels repeated


# pdf("CTD.pdf", height = 15, width=10)
# plots <- dlply(iS, ~variable, function(x)
#         {
#         ggplot(x, aes(x=distance, y=Depth.m)) +
#         #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
#         geom_raster(aes(fill=value, na.rm=T))+
#         stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
#         scale_fill_gradientn(paste(x$variable), colours=spectral(), na.value=NA) +
#         scale_x_continuous("Distance from shore", expand=c(0,0)) +
#         scale_y_reverse("Depth (m)", expand=c(0,0)) +
#         opts
#         })
#
#
# do.call(grid.arrange, c(plots,list(ncol=1)))
# dev.off()


# longer but better way
data <- iS[which(iS$variable == "Salinity.PPT"), ]
p1 <- ggplot(data, aes(x=distance, y=Depth.m)) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(data$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("", expand=c(0,0)) +
        scale_y_reverse("", expand=c(0,0)) +
        opts

data <- iS[which(iS$variable == "Temp.celsius"), ]
p2 <- ggplot(data, aes(x=distance, y=Depth.m)) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(data$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("", expand=c(0,0)) +
        scale_y_reverse("Depth (m)", expand=c(0,0)) +
        opts

data <- iS[which(iS$variable == "Fluoro.volts"), ]
p3 <- ggplot(data, aes(x=distance, y=Depth.m)) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(data$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("", expand=c(0,0)) +
        scale_y_reverse("", expand=c(0,0)) +
        opts

data <- iS[which(iS$variable == "Oxygen.ml.l"), ]
p4 <- ggplot(data, aes(x=distance, y=Depth.m)) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(data$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("Distance from shore", expand=c(0,0)) +
        scale_y_reverse("", expand=c(0,0)) +
        opts

pdf("ctd-clean.pdf", height = 13, width=9)
grid.arrange(p1, p2, p3, p4, ncol=1)
dev.off()


 }


#------------------------------------------------------
#                PLOT TAXA ABUNDANCES
#------------------------------------------------------


# select variables
colnames(biophy)
biophyplot <- biophy# [, c(1, 2, 4, 6:11, 13:16, 18:20, 23, 25:33, 47, 49, 93, 91)]

sal <- iS[which(iS$variable == "Salinity.PPT"), ]
fluo <- iS[which(iS$variable == "Fluoro.volts"), ]
temp <- iS[which(iS$variable == "Temp.celsius"), ]
oxy <- iS[which(iS$variable == "Oxygen.ml.l"), ]

# load pca data to plot
load("pca-physics.RData")


# RADIOLARIANS
# --------------------------------------------------------------------
 {
     
data <- biophyM[which(biophyM$taxa == "radiolarians_solitarian_dark"), ]
p_rd_val <- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("radiolarians dark m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "radiolarians_solitarian_other"), ]
p_rs_val <- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("radiolarians sol m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    

data <- biophyM[which(biophyM$taxa == "radiolarians_acantharia"), ]
p_aca_val <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("acantharians m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "radiolarians_colony"), ]
p_rc_val <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("radiolarians col m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    
pdf("plot/rads_validated.pdf", height = 15, width=10)
grid.arrange(p_rd_val, p_rs_val, p_aca_val, p_rc_val, ncol=1)
dev.off()


 }


# JELLYFISHES
# --------------------------------------------------------------------
 {

data <- biophyM[which(biophyM$taxa == "jelly_arctapod"), ]
p_jarc <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Arctap m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


# data <- biophyM[which(biophyM$taxa == "jelly_solmis"), ]
# p_jsol <- ggplot() +
#     geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
#     stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
#     geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
#     scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
#     #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
#     scale_size(expression(paste("Solmis m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
#     scale_x_continuous("", expand=c(0,0)) +
#     scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
#     opts


data <- biophyM[which(biophyM$taxa == "jelly_other"), ]
p_joth <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Others m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts

p_joth_pca <- ggplot() +
    geom_raster(data=pca, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=pca, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Others m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "ephyrae"), ]
p_eph <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Ephyrae m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


pdf("jellies_sal_validated.pdf", height = 15, width=10)
grid.arrange(p_jarc, p_eph, p_joth, p_jsol, ncol=1)
dev.off()

pdf("jellies_other_pca_validated.pdf", height = 9, width=10)
grid.arrange(p_joth, p_joth_pca, ncol=1)
dev.off()

 }


# OTHER GELATINOUS
# --------------------------------------------------------------------
 
 {

data <- biophyM[which(biophyM$taxa == "doliolids"), ]
p_dol_val <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Doliolids m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



p_dol_val_log <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=log(abund.m3+1)))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Doliolids m"^"-3")), limits = c(1, max(log(data$abund.m3+1))), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



pdf("doliolids_val_sal.pdf", height = 9, width=10)
grid.arrange(p_dol_val, p_dol_val_log, ncol=1)
dev.off()



# SIPHONOPHORES
# --------------------------------------------------------------------

data <- biophyM[which(str_detect(biophyM$taxa,  "siphos") == T), ]

p_sipho <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Calycophores m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



p_sipho_log <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=log(abund.m3+1)))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Calycophores m"^"-3")), limits = c(1, max(log(data$abund.m3+1))), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



pdf("sipho_log_sal.pdf", height = 9, width=10)
grid.arrange(p_sipho, p_sipho_log, ncol=1)
dev.off()




# Ctenophores
# --------------------------------------------------------------------


data <- biophyM[which(biophyM$taxa == "cteno_mertens"), ]

p_cteno <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Ctenophores m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts

pdf("cteno_sal.pdf", height = 6, width=10)
grid.arrange(p_cteno, ncol=1)
dev.off()



# APPENDICULARIANS
# --------------------------------------------------------------------

data <- biophyM[which(biophyM$taxa == "append_oikop"), ]
p_oikop <- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Oikop m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



data <- biophyM[which(biophyM$taxa == "append_fritill"), ]
p_frit <- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("fritill m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



pdf("append_fluo.pdf", height = 9, width=10)
grid.arrange(p_oikop, p_frit, ncol=1)
dev.off()


 }


# CRUSTACEANS
# --------------------------------------------------------------------

  {

# COPEPODS
# --------------------------------------------------------------------

  

data <- biophyM[which(biophyM$taxa == "crust_copepods"), ]
p_cop_val <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Copepods m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "crust_copepods_calanus"), ]
p_calanus <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Calanus m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "crust_copepods_euchaeta"), ]
p_euch <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("euchaeta m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "crust_copepods_other"), ]
p_copOth <- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Others m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


pdf("cop.pdf", height = 15, width=10)
grid.arrange(p_cop, p_calanus, p_euch, p_copOth, ncol=1)
dev.off()





# CUSTACEANS OTHERS
# --------------------------------------------------------------------

data <- biophyM[which(biophyM$taxa == "crust_amphipods"), ]
p_amph <- ggplot() +
    geom_raster(data=temp, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=temp, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Amphipods m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    
    
data <- biophyM[which(biophyM$taxa == "crust_larvae"), ]
p_larv <- ggplot() +
    geom_raster(data=pca, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=pca, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Crust larv m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts



data <- biophyM[which(biophyM$taxa == "polychaetes"), ]
p_poly <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Polychaetes m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts

pdf("others.pdf", height = 12, width=10)
grid.arrange(p_amph, p_larv, p_poly, ncol=1)
dev.off()


  }


# OTHER GROUPS
# --------------------------------------------------------------------

 {

# Chaetognaths
# --------------------------------------------------------------------


data <- biophyM[which(biophyM$taxa == "chaetognaths"), ]
p_chaeto <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Chaetognaths m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts

pdf("chaeto_sal.pdf", height = 6, width=10)
grid.arrange(p_chaeto, ncol=1)
dev.off()



# Shrimps
# --------------------------------------------------------------------

data <- biophyM[which(biophyM$taxa == "shrimps"), ]
p_shrimps <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Shrimps m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts

pdf("shrimps_sal.pdf", height = 6, width=10)
grid.arrange(p_shrimps, ncol=1)
dev.off()


# PTEROPODS
# --------------------------------------------------------------------


data <- biophyM[which(biophyM$taxa == "pteropods_cavol"), ]
p_cavol <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Cavolina m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts

data <- biophyM[which(biophyM$taxa == "pteropods_creseis"), ]
p_cres <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Creseis m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


data <- biophyM[which(biophyM$taxa == "pteropods_other"), ]
p_pteOth <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Other m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


pdf("pte_sal.pdf", height = 12, width=10)
grid.arrange(p_cavol, p_cres, p_pteOth, ncol=1)
dev.off()



# PHYTO
# --------------------------------------------------------------------

data <- biophyM[which(biophyM$taxa == "phyto_diatom_chains"), ]
p_diat <- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Diatoms chains m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    
    
data <- biophyM[which(biophyM$taxa == "phyto_trichod"), ]
p_trich<- ggplot() +
    geom_raster(data=fluo, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=fluo, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Trichodesm m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    

pdf("phyto_fluo.pdf", height = 9, width=10)
grid.arrange(p_diat, p_trich, ncol=1)
dev.off()


 }    

    
# FISH AND FISH-LIKE
# --------------------------------------------------------------------

 {
    
data <- biophyM[which(biophyM$taxa == "fish_larvae"), ]
p_fish<- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Fish larvae m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    
data <- biophyM[which(biophyM$taxa == "fish_like_high_conf"), ]
p_fish_hc <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Fish-hConf m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts
    
    
data <- biophyM[which(biophyM$taxa == "fish_like"), ]
p_fishlike <- ggplot() +
    geom_raster(data=sal, aes(x=distance, y=-Depth.m, fill=value, na.rm=T))+
    stat_contour(data=sal, aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.8, bins=5, size=0.4, na.rm=TRUE) +
    geom_point(data=data, aes(x=distanceFromVlfr, y=-DepthBin, size=abund.m3))+
    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
    #scale_size(bquote("Ind m"^"-3"), limits = c(1, max(x$value)), range = c(1, 12))+
    scale_size(expression(paste("Fish-like m"^"-3")), limits = c(1, max(data$abund.m3)), range = c(1, 12))+
    scale_x_continuous("", expand=c(0,0)) +
    scale_y_continuous("", limits = c(-103, 5), expand=c(0,0)) +
    opts


pdf("fish_sal_all.pdf", height = 12, width=10)
grid.arrange(p_fish, p_fish_hc, p_fishlike, ncol=1)
dev.off()


 }


#---------------------------------------------
#           PLOT SHIP TRAJECTORY 
#---------------------------------------------
 {


# YO-YOs

# overplot ISIIS trajectory over physical data
# Get physical data before interpolation
updw <- phy

# select the processed casts
processed <- c(2, 7, 9, 17, 25, 33, 49)
updw[which(updw$cast %in% processed), "Done"] <- "Processed"
updw[which(is.na(updw$Done)), "Done"] <-  "To be processed"


pdf("yoyos.pdf", width = 9, height = 4)
ggplot() +
        geom_raster(aes(x=distance, y=Depth.m, fill=value), data= sal, na.rm=T, ) +
        stat_contour(aes(x=distance, y=Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE, data=sal) +
        geom_line(aes(x=distanceFromVlfr, y=Depth.m, linetype=Done, group=cast), size=0.5, data=updw) +
        scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
        scale_x_continuous("Distance from shore (nm)", expand=c(0,0)) +
        scale_y_reverse("Depth (m)", expand=c(0,0)) +
        scale_linetype_manual("", values=c(1, 3)) + 
        opts
dev.off()




# read coastline to plot the trajectories and check
coast <- read.csv(str_c(dir, "map/gshhg_coteazur_i.csv"))
load(str_c(dir, "map/coast_bathy.RData"))

# read stations position
station <- read.csv(str_c(dir, "/nets/station-regent-visufront.csv"), header=T, sep=";")
head(station)

# compute lat and lon for stations
latBits <- str_split_fixed(station$lat_out, fixed("."), 2)
station$lat <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
station$lat <- 43 + station$lat/60
lonBits <- str_split_fixed(station$lon_out, fixed("."), 2)
station$lon <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
station$lon <- 7 + station$lon/60

# plot station position
ggplot(mapping=aes(x=lon, y=lat)) + 
	geom_polygon(data=coast, fill="grey60") + 
	geom_point(aes(color=date), data=station, size=3) + 
	geom_text(data=station, label=station$station_nb, vjust=2) + 
	geom_point(data=Boussole, color="blue") + 
	coord_map(xlim=c(7,8.1), ylim=c(43.2,43.75))



# plot boat traj + stations (from drifter-plot.R)
# read ship trajectory from ts
filenames <- list.files(str_c(dir, "/TS/"))

s <- adply(filenames, 1, function(x) {
	s <- read.ts(str_c(dir, "/TS/",x))
	return(s)
	}, .progress="text")


# Plot ship trajectory and stations
ggplot(mapping=aes(x=lon, y=lat)) + geom_contour(aes(z=-z, x=x, y=y), colour="gray80", data=bathyDF, size=0.3) + geom_polygon(fill="gray25", data=coast, aes(x=lon, y=lat)) + geom_point(data=station, size=4, colour= "gray40") + geom_path(size=0.35, na.rm=T, data=s, colour="gray20") + scale_x_continuous("Longitude", expand=c(0,0)) + scale_y_continuous("Latitude", expand=c(0,0)) +scale_color_discrete("") + coord_quickmap(xlim=c(6.9, 8.05), ylim=c(43.24, 43.75)) +theme_bw() + opts



# Add glider data that to complete the profile
#----------------------------------------------
g <- read.csv("~/Desktop/PhD/VISUFRONT/ISIIS/glider/glier-d25.csv", sep=",", header=T)
head(g)

min(g$distance)
g <- g[, c("ilon", "ilat", "distance")]
g <- rename(g, c('ilon'='lon', 'ilat'='lat'))

s24_25 <- s[which(s$dateTime > as.POSIXct("2013-07-24 15:00:00") & s$dateTime < as.POSIXct("2013-07-25 05:00:00")), ]

ggplot(mapping=aes(x=lon, y=lat)) + geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) + geom_polygon(fill="gray25", data=coast, aes(x=lon, y=lat)) + geom_path(size=0.4, na.rm=T, data=s24_25) + geom_point(size=2, data=g)+ scale_x_continuous("Longitude", expand=c(0,0)) + scale_y_continuous("Latitude", expand=c(0,0)) + coord_quickmap(xlim=c(6.8, 8.1), ylim=c(43.2, 43.75)) + theme_bw()




 }


#---------------------------------------------------------
#          COMPUTE ABUNDANCE PER VOL PLKT NETS
#---------------------------------------------------------
 {

# Read abundance at stations
fl_nets <- read.csv(str_c(dir, "/nets/visufront-fish-larvae.csv"), sep=";", h=T)
head(fl_nets)

# COUNTS
sum(fl_nets$larvae_nb) # 675 fish larvae
sum(station$larvae) # 677 -> 2 cephalopodidae (?)
sum(station$eggs) # 548 eggs
length(unique(fl_nets$order)) # 13 orders
length(unique(fl_nets$family)) # 29 families

# species nb
species <- unique(str_c(str_c(fl_nets$genera, fl_nets$sp, sep=" "), str_c(" (", fl_nets$descriptor, ")")))
length(species)  # 47 species

# Read species list with their related habitats
#write.table(species, "species-only.csv", row.names=F, col.names=F)
species <- read.csv(str_c(dir, "nets/species.csv"), h=T, sep=";")
species$species <- as.character(species$species)

# add a column with the species name in the fl_nets
fl_nets$species <- str_c(str_c(fl_nets$genera, fl_nets$sp, sep=" "), str_c(" (", fl_nets$descriptor, ")"))

# add a column to the stations to join it to the fish abundance
station$station <- str_c("station_", station$station_nb)
d <- join(fl_nets, station, by="station")
d <- join(d, species)


#------------------------------------------
# Compute abundance per volume unit
#------------------------------------------

# 1st method = standardized abundance by volume using volumeter. 
#-------------------------------------------------------------------
# NB : Volume unit remains unknown

# compute abundance stadardized per volume unit for stations from day 24 only
# NB : volume unit is unkown here, only based on the volumeter diff between in and out, but it already makes absolute data comparable
abund.vol <- ddply(d, .(station, habitat, date), function(x){
    data.frame(sum=sum(x$larvae_nb), 
    abundance = sum(x$larvae_nb) / (x$vol_out[1]-x$vol_in[1]), 
    order = length(unique(x$order)), 
    family = length(unique(x$family)), 
    species = length(unique(x$sp)), 
    lat=x$lat[1], lon=x$lon[1], date=x$date[1])
    })

head(abund.vol)



# 2nd method = compute sampled volume using speed and distance
#------------------------------------------------------------------

# select variables  of lon and lat
station <- ddply(d, ~station, function(x){
    data.frame(lon=unique(x$lon), lat=unique(x$lat),
    latin=unique(x$lat_in), lonin=unique(x$lon_in),
    latout=unique(x$lat_out), lonout=unique(x$lon_out), 
    time_in = str_c("2013-07-18 ", unique(x$time_start)),
    time_out = str_c("2013-07-18 ", unique(x$time_end)),
    larvae_tot = unique(x$larvae)) })


# compute and format lat and lon for stations
latBits <- str_split_fixed(station$latout, fixed("."), 2)
station$lat_out <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
station$lat_out <- 43 + station$lat_out/60
latBits <- str_split_fixed(station$latin, fixed("."), 2)
station$lat_in <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
station$lat_in <- 43 + station$lat_in/60

lonBits <- str_split_fixed(station$lonout, fixed("."), 2)
station$lon_out <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
station$lon_out <- 7 + station$lon_out/60
lonBits <- str_split_fixed(station$lonin, fixed("."), 2)
station$lon_in <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
station$lon_in <- 7 + station$lon_in/60


# compute distance and time per station
dist_time <- ddply(station, ~station, function(x){
    difftime <- difftime(as.POSIXct(x$time_out), as.POSIXct(x$time_in))
    dist.m <- geodDist(lat1=x$lat_in, lon1=x$lon_in, lat2=x$lat_out, lon2=x$lon_out, alongPath=F)
    time <- as.numeric(difftime) * 60
    speed <- dist.m / time * 1000
    vol.s <- speed * pi * 1^2  # pi.r^2 = 0.78 m^2
    vol.m3 <- vol.s * time 
    return(data.frame(vol.m3))})


# Very bad but some positions may be wrong so I fill replace them by possibly consistent values
dist_time[which(dist_time$station =="station_7"), "vol.m3"] <- 350
dist_time[which(dist_time$station =="station_2"), "vol.m3"] <- 330
dist_time[which(dist_time$station =="station_15"), "vol.m3"] <- 650


# join station sampled volume and larval abundances
abund.m3 <- join(station, dist_time)
abund.m3

abund.m3$abund.m3 <- abund.m3$larvae_tot / abund.m3$vol.m3
abund.m3$abund.m3.norm <- abund.m3$abund.m3 / max(abund.m3$abund.m3)
abund.m3$data <- "Plankton nets"


 }



# --------------------------------------------------------------------
#       COMPARE ISIIS AND PLANKTON NETS ABUNDANCES
# --------------------------------------------------------------------
 {

# get abundance from ISIIS data
lf_isiis <- biophy[, names(biophy) %in% c("fish", "lat", "lon", "cast")]

# integrate larval fish abundance over profiles
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# x$fish is abund.m3 and not number of larvae
lf_isiisSummary <- ddply(lf_isiis, ~cast, function(x) {
data.frame(sum = sum(x$fish), lat = x$lat[1], lon = x$lon[1], data="ISIIS")
})


# compute integrated volume per cast
volisiis <- ddply(vol, ~cast, function(x){
    data.frame(vol = sum(unique(x$vol.m3)))
})
    
# join integrated abundance and integrated volume
isiisvol <- join(lf_isiisSummary, volisiis)
isiisvol$abund.m3 <- isiisvol$sum / isiisvol$vol
isiisvol$abund.m3.norm <- isiisvol$abund.m3 / max(isiisvol$abund.m3)
    
isiisvol  <- rename(isiisvol, c("sum"="larvae.tot", "cast"="station"))
    
    
    
# select the station we want
abund.m3 <- abund.m3[-which(abund.m3$station %in% c("station_19", "station_20")),]


# RBIND ISIIS AND PLANKTON NET SAMPLES
fishab <- rbind(isiisvol, ab) 
    
p <- ggplot() + 
#	geom_raster(aes(fill=-z, x=x, y=y), data=bathyDF) + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray80", data=bathyDF, size=0.3) + 
	geom_polygon(fill="gray25", data=coast, aes(x=lon, y=lat)) +
    geom_point(aes(x=lon, y=lat, size=abund.norm, colour=data), data=fishab) +
    scale_x_continuous("Longitude", expand=c(0,0)) + 
	scale_y_continuous("Latitude", expand=c(0,0)) +
    scale_size("Relative abundance", range=c(1,10))+#, breaks = c(1, 5, 10)) +
	coord_quickmap(xlim=c(6.9, 8.05), ylim=c(43.24, 43.75)) +
	theme_bw() #+ opts

pdf("isiis-nets-comp.pdf", width=9 ,height=6)
p
dev.off()





# !!!!!!!!!!!! 
# ADD A PLOT WITH THE 2 ABUNDANCES OF PLANKTON NETS CALCULATED (M3 and VOL)









# Plots of larval fish STANDARDIZED ABUNDANCE
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=sum), colour="gray20", data=abund.vol[which(abund.vol$date %in% c("2013-07-18", "2013-07-19", "2013-07-26")), ]) + 
	scale_size_area("Relative Abund") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))



# For larval fish diversity NUMBER OF SPECIES
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=species), colour="gray20", data=abund.vol) + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Species", breaks=c(2, 5, 10)) +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))



# For larval fish diversity OTHER THAT COASTAL (pelagic, mesopelagic and benthopelagic)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=abundance), colour="gray20", data=abund[-which(abund$habitat %in% "coastal"), ]) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Abundance offshore") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))


# For larval fish diversity OTHER THAT COASTAL (pelagic, mesopelagic and benthopelagic)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=abundance), colour="gray20", data=abund[which(abund$habitat %in% "mesopelagic"), ]) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Abundance Meso") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))




# For larval fish diversity OTHER THAT COASTAL (pelagic, mesopelagic and benthopelagic)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=order), colour="gray20", data=abund[which(abund$habitat %in% "coastal"), ]) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Abundance coastal") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))

dev.off()		


 }



#---------------------------------------------
#                   TODO
#---------------------------------------------

# ADD A PLOT WITH THE 2 ABUNDANCES OF PLANKTON NETS CALCULATED (M3 and VOL)

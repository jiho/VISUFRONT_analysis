# -----------------------------------------------------------
#
# Process ISIIS data for the Ocean Science Meeting 2014
# Robin Faillettaz, 2014-02-14
#
# -----------------------------------------------------------
#
# Working on transect cross - front 04
# It's a night transect
# Profil validated : 01, 03, 05, 09, 13, 17, 25 (?)
#
# A .Rprofile file should be set personnaly to the directory where the data are stored (dropbox)

# For now locate data from the dropbox repo manually
dir <- "/Users/faillettaz/Dropbox/robin/visufront-data"
dir <- "/Users/faillettaz/Dropbox/visufront-data/"


library("plyr")
library("ggplot2")
library("stringr")
library("reshape2")
library("grid")
library("gridExtra")
library("akima")
library("fields")


source("data/lib_zooprocess.R")
source("lib_zooprocess.R")
source("data/lib_plot.R")


files <- list.files(dir, full = T)
files <- list.files(str_c(dir, "zooprocess/"), full = T)
#pid <- read.pid(files[2])

dat1 <- files[which(str_detect(files, "_dat1.txt") == T)]
datfiles <- files[which(str_detect(files, "_datfile.txt") == T)]

pids <- NA
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
rbind(pids, pid)

}, .progress="text")  
pids <- pids[-1, -1]

# check for NAs
pids[which(is.na(pids)), ]

# delete them if any
if (dim(pids[which(is.na(pids)), ])[1] > 0) {
   pids <- pids[-which(is.na(pids)), ] 
}






#--------------------------------------------
#       CLEAN NAMES AND MERGE GROUPS
#--------------------------------------------


# Check "not_found" category
pids[which(pids$Valid == 'not_found'), ]

# correspond to images lost during sorting, get rid of them
if (length(pids[which(pids$Valid == 'not_found'), "Valid"] > 0)) {
    pids <- pids[-which(pids$Valid == 'not_found'), ]
}


# Merge SIPHOS
pids[which(pids$Valid %in% c("sipho_tail", "sipho_round")), "Valid"] <- "sipho"
unique(pids$Valid)


# Merge EPHYRAE
pids[which(pids$Valid %in% c("ephyrae_side")), "Valid"] <- "ephyrae"
sort(unique(pids$Valid))


# Merge polychaetes
pids[which(pids$Valid %in% c("polychaets")), "Valid"] <- "polychaetes"
sort(unique(pids$Valid))


# Merge shrimps
pids[which(pids$Valid %in% c("shrimp_large")), "Valid"] <- "shrimps"
sort(unique(pids$Valid))

# change radiolarians names
pids[which(pids$Valid %in% c("radiolarians")), "Valid"] <- "radiolarian_sol"
sort(unique(pids$Valid))
pids[which(pids$Valid %in% c("radiolarian_col_rings")), "Valid"] <- "radiolarian_rings"
sort(unique(pids$Valid))
pids[which(pids$Valid %in% c("radiolarian_colony")), "Valid"] <- "radiolarian_col"
sort(unique(pids$Valid))
pids[which(pids$Valid %in% c("radiolarians_dark")), "Valid"] <- "radiolarian_dark"
sort(unique(pids$Valid))



# Merge groups
pids$groups <- pids$Valid 
# change radiolarians names
pids[which(str_detect(pids$Valid, "radiolarian") == T), "groups"] <- "radiolarians"
pids[which(str_detect(pids$Valid, "fish") == T), "groups"] <- "fish_like"
pids[which(pids$Valid %in% c("crust_larvae", "crustaceans", "copepods")), "groups"] <- "small_crustaceans"

sort(unique(pids$groups))






# set depth BIN 
binsize <- 0.5
pids$DepthBin <- round_any(pids$Depth, binsize)

# get raw abundance
bio <- ddply(pids, ~Valid+Label+DepthBin, function(x) {sum(na.omit(x$Valid==paste(x$Valid)))})
bio <- rename(bio, c("V1" = "Abund"))

# check for NAs
bio[which(is.na(bio$Valid)), ]

if (dim(bio[which(is.na(bio)), ])[1] > 0) {
    bio <- bio[-which(is.na(bio)), ]
}




# join with depth bin to have a line for each depth
# Set all possibilities of depth, valid and label and add corresponding abundances
grid <- expand.grid(unique(bio$Valid), unique(bio$Label), unique(bio$DepthBin))
head(grid)
names(grid) <- c("Valid", "Label", "DepthBin")

# join the grid with abundances
bioFull <- join(grid, bio)


# replace NAs by 0
bioFull[which(is.na(bioFull$Abund)), "Abund"] <- 0 
length(which(is.na(bioFull)))  # if 0 --> OK


# get number of a certain group
sum(bioFull[which(bioFull$Valid == "aggregates"), "Abund"])


# get cast number from the profil name (vignette label)
bioFull$cast <- as.numeric(str_split_fixed(bioFull$Label, fixed("_"), 3)[, 3]) * 2 - 1


# plot everything
ggplot(bioFull[which(bioFull$Valid %in% c("fish", "fish_like", "chaetognaths", "sipho")), ]) + 
geom_point(aes(x = Abund, y = -DepthBin, colour = Label)) + 
geom_path(aes(x = Abund, y = -DepthBin, colour = Label)) + 
facet_grid(.~Valid, scales="free_x")




#--------------------------------------------
#       Get sampled volumn per bin
#--------------------------------------------

vol <- adply(datfiles, 1, function(x) {
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
    img$cast <- as.numeric(str_split_fixed(x, fixed("_"), 4)[, 3]) * 2 -1
    
    
    return(img)
    
    }, .progress="text")
    vol <- vol[, -1]
    
head(vol)



#--------------------------------------------
#       Merge volumns with abundance
#--------------------------------------------

head(bioFull)
head(vol)
bioFull <- join(bioFull, vol)
bioFull$abund.m3 <- bioFull$Abund / bioFull$vol.m3



#---------------------------------------------------------
#         Process Physical Data
#---------------------------------------------------------

# Read physical data from transect 4
phy <- read.csv("transects/cross_current_4/isiis.csv", header=T, sep=",")
phy <- read.csv(str_c(dir, "/isiis.csv"), header=T, sep=",")
head(phy)
# delete first line (data from previous transect)
phy <- phy[-1, ]

# select only data from upcasts --> Keep first down cast below 25m to improve interpolation
d <- phy[which(phy$down.up %in% "up" | phy$down.up %in% "down" & phy$Depth.m > 28), ]


# check if seems ok
ggplot(d, aes(x=distanceFromVlfr, y=-Depth.m, colour=Salinity.PPT))  + geom_point() # Yes !




#-----------------------------------
#        GLIDER DATA
#-----------------------------------

# Add glider data that to complete the profile
g <- read.csv("~/Desktop/PhD/VISUFRONT/ISIIS/glider/glier-d25.csv", sep=",", header=T)
head(g)


# check path
ggplot(g) + geom_path(aes(x=distance, y=ipressure, colour= substring(dateTime, 1, 10)))

# get data from cast 1 to 4
g[which(round(g$distance) == 7), ]
min(g$distance)


# select data to keep 
# keep only variable of interest
g <- g[which(g$ipressure < 103 & g$distance < 7.2), names(g) %in% c("distance", "ipressure", "salinity")]


d_short <- d[ , names(d) %in% c("distanceFromVlfr", "Depth.m", "Salinity.PPT")]
head(d_short)

colnames(g) <- names(d_short)
head(g)

dg <- rbind(d_short, g)
head(dg)


# compute interpolation 
# interpolate all variables
dm <- melt(dg, id.vars=c("Depth.m", "distanceFromVlfr"), measure.vars=c("Salinity.PPT"))#, "Temp.C", "Fluoro.volts", "Oxygen.ml.l"))


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



#-----------------------------------
#        IF NO GLIDER DATA 
#-----------------------------------

# compute interpolation 
# interpolate all variables
dm <- melt(d, id.vars=c("Depth.m", "distanceFromVlfr"), measure.vars=c("Salinity.PPT"))#, "Temp.C", "Fluoro.volts", "Oxygen.ml.l"))



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

# delete NAs from the previous interpolation
i2 <- na.exclude(di)

# Select variables
x <- i2$distance
y <- i2$Depth.m
z <- i2$value




# Second interpolation for all variables
di2 <- ddply(di, ~variable, function(x) {
    x <- na.omit(x)
    xi <- interp.smooth(x=x$distance, y=x$Depth.m, z=x$value, x.step = 0.1, y.step = 0.1)
}, .progress="text")


di2 <- rename(di2, c("x"="distance", "y"="Depth.m"))

plots <- dlply(di2, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=-Depth.m)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_raster(aes(fill=value, na.rm=T)+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
        })

do.call(grid.arrange, c(plots,list(ncol=1)))



# overplot ISIIS trajectory over physical data
# Get physical data before interpolation
updw <- phy

# select the processed casts
processed <- c(1, 5, 9, 17, 25, 33, 49)
updw <- updw[which(updw$cast %in% processed), ]

# select a single variable to go faster
sal <- di2[which(di2$variable == "Salinity.PPT"), ]

ggplot() +
geom_raster(aes(x=distance, y=-Depth.m, fill=value), data= sal, na.rm=T, ) +
stat_contour(aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE, data=sal) +
geom_path(aes(x=distanceFromVlfr, y=-Depth.m, group=cast), size=0.6, data=updw) +
scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
scale_x_continuous(expand=c(0,0)) +
scale_y_continuous(expand=c(0,0))




#--------------------------------------------------
#           Merge physical and biological
#--------------------------------------------------



# NB : Bio data need to be localised in x, y, z --> I will bin the depth by profile
# and get the corresponding x and y


head(bioFull)


# NB : Bio data need to be localised in x, y, z --> I will bin the depth by profile
# and get the corresponding x and y

# select only physical data from the profil 05 (cast 9)
biophy <- ddply(bioFull, ~cast, function(x) {
    
    phyx <- phy[which(phy$cast == paste(x$cast)), ]
    # Bin depth = select middle of each meter
    # set depth BIN 
    bin <- seq(from=1, to = ceiling(max(phyx$Depth.m)), by=binsize)
    
    # get index of the min of each point --> NEED TO BE DONE BY PROFILE
    min <- adply(bin, 1, function(m){
        which.min(abs(phyx$Depth.m-m))
        })

    # Merge biology and physics
    biox <- bioFull[which(bioFull$cast == paste(x$cast)), ]


    # PROBLEMS WITH DCAST BUT ELSE IT'S SEEMS TO WORK

    dC <- dcast(biox, DepthBin ~ Valid, value.var = "abund.m3")
    length(which(is.na(dC)))

    # select only the physical data corresponding to the center of each bin
    phyxMin <- phyx[min$V1, ]
    phyxMin$DepthBin <- round_any(phyxMin$Depth.m, binsize)

    if (phyxMin$DepthBin[length(phyxMin$DepthBin)] == 102) {
        phyxMin$DepthBin[length(phyxMin$DepthBin)] <- 103
    }

    dim(phyxMin)
    dim(dC)

    # join the 2
    biophyx <- join(dC, phyxMin)
    
    return(biophyx)
    
}, .progress="text")

# NB :Warnings problably because some profiles don't go as deep as or deeper that 103m

head(biophy)
dim(biophy)
length(unique(biophy$cast))
table(biophy$DepthBin)


ggplot() + 
geom_raster(aes(x=distance, y=-Depth.m, fill=value), data= sal, na.rm=T, ) +
stat_contour(aes(x=distance, y=-Depth.m, z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE, data=sal) +
#geom_path(aes(x=distanceFromVlfr, y=-Depth.m, group=cast), size=0.6, data=phy) +
geom_point(data=biophy[-which(biophy$fish==0), ], aes(x=distanceFromVlfr, y=-DepthBin, size=fish), group=cast)+
scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
scale_x_continuous(expand=c(0,0)) +
scale_y_continuous(expand=c(0,0))




#------------------------------------------------------
#         Get physical data from distribution 
#------------------------------------------------------

head(biophy)
xy <- biophy[, c("distanceFromVlfr", "Depth.m")]

# interpolate
i <- interp(x=x, y=y, z=z, xo=seq(0, max(x), by=0.1), yo=seq(0, max(y), by=0.1),  linear = T, duplicate = "mean", extrap=F)


grid <- data.frame(x=seq(0, max(x), by=0.1), y=seq(0, max(y), length=length(seq(0, max(x), by=0.1))))

smooth <- image.smooth(i, grid=grid, theta=0.28)
out <- melt(smooth$z, varnames=c("x","y"))
out$x <- smooth$x[out$x]
out$y <- smooth$y[out$y]

out <- rename(out, c("x"="distance", "y"="Depth.m"))


ggplot(out, aes(x=distance, y=-Depth.m))+
#geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
geom_raster(aes(fill=value, na.rm=T))+
stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
geom_point(data=biophy[-which(biophy$fish==0), ], aes(x=distanceFromVlfr, y=-DepthBin, size=fish), group=cast)+
scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
scale_size(expression(paste("Fish larvae.m"^"-3")))+
scale_x_continuous("Distance from shore", expand=c(0,0)) +
scale_y_continuous("Depth (m)", expand=c(0,0))




plots <- dlply(di2, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=-Depth.m)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
        })

do.call(grid.arrange, c(plots,list(ncol=1)))




biophy$sal <- interp.surface(smooth, xy)

biophy[, c("sal", 'Salinity.PPT', 'Depth.m')]

#---------------------------------------------
#                   TODO
#---------------------------------------------

Compute total larval abundance per station to compare with plankton nets
interp.surface(xy, xyz)





#---------------------------------------------
#           Plot ship trajectory 
#---------------------------------------------

# read coastline to plot the trajectories and check
coast <- read.csv("map/cote_azur.csv")
load("map/coast_bathy.RData")
load(str_c(dir, "/map/coast_bathy.RData"))

# read stations position
station <- read.csv("Plankton-nets/station-regent-visufront.csv", header=T, sep=";")
station <- read.csv(str_c(dir, "/nets/station-regent-visufront.csv"), header=T, sep=";")
head(station)

# compute lat and lon for stations
latBits <- str_split_fixed(station$lat_out, fixed("."), 2)
station$lat <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
station$lat <- 43 + station$lat/60
lonBits <- str_split_fixed(station$lon_out, fixed("."), 2)
station$lon <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
station$lon <- 7 + station$lon/60

# Add Boussole as the reference point
Boussole <- data.frame(lat=43.38, lon=7.83)

p <- ggplot(mapping=aes(x=lon, y=lat)) + 
	geom_polygon(data=coast, fill="grey60") + 
	geom_point(aes(color=date), data=station, size=3) + 
	geom_text(data=station, label=station$station_nb, vjust=2) + 
	geom_point(data=Boussole, color="blue") + 
	coord_map(xlim=c(7,8.1), ylim=c(43.2,43.75))
print(p)



# plot boat traj + stations (from drifter-plot.R)
# read ship trajectory from ts
filenames <- list.files("TS/")
filenames <- list.files(str_c(dir, "/TS/"))

source("lib_process.R")

s <- adply(filenames, 1, function(x) {
	s <- read.ts(str_c("TS/",x))
	s <- read.ts(str_c(dir, "/TS/",x))
	return(s)
	})

ggplot(mapping=aes(x=lon, y=lat)) + 
	#geom_raster(aes(fill=-z, x=x, y=y), data=bathyDF) + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) + 
	geom_polygon(fill="gray25", data=coast, aes(x=lon, y=lat)) +
	geom_path(size=0.4, na.rm=T, data=s) + # ship track 
	geom_point(aes(color=date), data=station, size=3) + 
	geom_text(data=station, label=station$station_nb, vjust=2) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="blue") +
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="blue", vjust=2, hjust=1, size=4) +
	scale_x_continuous("Longitude", expand=c(0,0)) + 
	scale_y_continuous("Latitude", expand=c(0,0)) +
	#scale_fill_gradient(name="Depth") +
	scale_color_discrete("") + 
	coord_quickmap(xlim=c(6.8, 8.1), ylim=c(43.2, 43.75)) +
	theme_bw()


    # Add glider data that to complete the profile
g <- read.csv("~/Desktop/PhD/VISUFRONT/ISIIS/glider/glier-d25.csv", sep=",", header=T)
head(g)

min(g$distance)
g <- g[, c("ilon", "ilat", "distance")]
g <- rename(g, c('ilon'='lon', 'ilat'='lat'))

s24_25 <- s[which(s$dateTime > as.POSIXct("2013-07-24 15:00:00") & s$dateTime < as.POSIXct("2013-07-25 05:00:00")), ]

    ggplot(mapping=aes(x=lon, y=lat)) + 
    	#geom_raster(aes(fill=-z, x=x, y=y), data=bathyDF) + 
    	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) + 
    	geom_polygon(fill="gray25", data=coast, aes(x=lon, y=lat)) +
    	geom_path(size=0.4, na.rm=T, data=s24_25) + # ship track 
        geom_point(size=2, data=g)+
    	scale_x_continuous("Longitude", expand=c(0,0)) + 
    	scale_y_continuous("Latitude", expand=c(0,0)) +
    	coord_quickmap(xlim=c(6.8, 8.1), ylim=c(43.2, 43.75)) +
    	theme_bw()


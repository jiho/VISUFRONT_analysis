# -----------------------------------------------------------
#
#   Spatial PCA  on physical data to localize the front
#              Robin Faillettaz, 2014-03-11
#
# -----------------------------------------------------------

# Working on transect cross - front 04
# It's a night transect

# For now locate data from the dropbox repo manually
dir <- "/Users/faillettaz/Dropbox/visufront-data/"


library("plyr")
library("ggplot2")
library("stringr")
library("reshape2")
library("grid")
library("gridExtra")
library("akima")
library("fields")
library("oce")   
library("FactoMineR")


source("lib_zooprocess.R")
source("lib_plot.R")
source("lib_process.R")



# set options for ploting
opts <- theme(axis.title.y= element_text(angle=90, vjust=0.5, size=18), axis.title.x= element_text(angle=0, vjust=0.5, size=18)) + 
theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=15), axis.text.y  = element_text(angle=0, vjust=0.5, size=15))


# Read physical data from transect 4
phy <- read.csv(str_c(dir, "ISIIShydro/transects/cross_current_4/isiis.csv"), header=T, sep=",")
head(phy)
# delete first line (data from previous transect)
phy <- phy[-1, ]

# Change names to same number of character for each (for plots)
phy <- rename(phy, c("Temp.C" = "Temp.celsius"))

# select only data from upcasts --> Keep first down cast below 25m to improve interpolation
phy25mUp <- phy[which(phy$down.up %in% "up" | phy$down.up %in% "down" & phy$Depth.m > 29), ]


# check if seems ok
ggplot(phy25mUp, aes(x=distanceFromVlfr, y=-Depth.m, colour=Salinity.PPT))  + geom_point() # Yes !


# compute interpolation 
# interpolate all variables
dm <- melt(phy25mUp, id.vars=c("Depth.m", "distanceFromVlfr"), measure.vars=c("Salinity.PPT", "Temp.celsius", "Fluoro.volts", "Oxygen.ml.l"))



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





#-------------------------------------------------------------
#       COMPUTE SPATIAL PCA ON PHYSICAL DATA
#-------------------------------------------------------------


#   First for raw data 
#-----------------------------

# Prepare data to compute PCA 
d <- dcast(di, Depth.m + distance ~ variable, value.var = "value")
ACP <- na.exclude(d)
ACP <- ACP[order(ACP$distance),]

# Compute PCA
res <- PCA(X=ACP[,-c(1,2)], graph=T)


# plot it per axis
ACP$pca <- res$ind$coord[,2]
ggplot(data=ACP, aes(x=distance, y=-Depth.m)) + geom_raster(aes(fill=pca)) + geom_contour(aes(z=pca), colour="white") + scale_fill_gradientn("PCA", colours = spectral()) + scale_x_continuous(name="Distance", expand=c(0,0)) + scale_y_continuous(name="Depth", expand=c(0,0))


#      Now on anomalies
#---------------------------------

# Compute mean profile for each variable
head(d)
dAvg <- ddply(na.omit(d), ~Depth.m, function(x){
    data.frame(
        sal = mean(x$Salinity.PPT), 
        temp = mean(x$Temp.celsius),
        fluo = mean(x$Fluoro.volts),
        oxy = mean(x$Oxygen.ml.l))
})


# Get anomalies of each variables along the transect per vertical profile
dN <- na.exclude(ddply(d[which(d$distance < 28), ], ~distance, function(x){
    data.frame(
        depth = x$Depth.m[-1],
        sal = x$Salinity.PPT[-1] - dAvg$sal,
        temp = x$Temp.celsius[-1] - dAvg$temp, 
        fluo = x$Fluoro.volts[-1] - dAvg$fluo,
        oxy = x$Oxygen.ml.l[-1] - dAvg$oxy)
}))


#   Check if ok 
#----------------------------
dM <- melt(dN, id.vars = c("depth", "distance"))
plots <- dlply(dM, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=depth)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(x$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("Distance from shore", expand=c(0,0)) +
        scale_y_reverse("Depth (m)", expand=c(0,0)) +
        opts
        })

do.call(grid.arrange, c(plots,list(ncol=1)))
# it's working


# compute the PCA with variables anomalies
pca <- dN[order(dN$distance),]
res <- PCA(X=pca[,-c(1,2)], graph=F)


# Set a factor variable to identified the quality of the pca to differentiate the areas (coastal, frontal and central)
pca$position <- NA
pca$position[which(pca$distance < 14)] <- "coastal"
pca$position[which(pca$distance > 21)] <- "central"
pca$position[which(is.na(pca$position))] <- "frontal"


# call the acp_plot function to plot the pca
plot_acp(x = res, colour = pca$position)


# save all different possibilites
pca$pca1 <- res$ind$coord[,1] # ~41%
pca$pca2 <- res$ind$coord[,2] # ~27%
pca$pca3 <- res$ind$coord[,3] # ~11%
pca$pca4 <- res$ind$coord[,4] # ~9%


# plot the output of the pca for each axis
allPCA <- melt(pca[, !names(pca) %in% c("sal", "temp")], id.vars = c("depth", "distance", "position"), measure.vars = grep("pca", names(pca), value=TRUE))
plots <- dlply(allPCA, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=depth)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_raster(aes(fill=value, na.rm=T))+
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(paste(x$variable), colours=spectral(), na.value=NA) +
        scale_x_continuous("Distance from shore", expand=c(0,0)) +
        scale_y_reverse("Depth (m)", expand=c(0,0)) +
        opts
        })

do.call(grid.arrange, c(plots,list(ncol=2)))




# Second interpolation over a finer grid for smoothing
#-----------------------------------------------------------


#   For axis 1 only
#-----------------------------

# create a new df and interpolate smooth from it
i2 <- data.frame(distance = dN$distance, Depth.m = dN$depth, value = pca$pca1)
di2 <- interp.smooth(x=i2$distance, y=i2$Depth.m, z=i2$value, x.step = 0.1, y.step = 0.1)
di2 <- rename(di2, c("x"="distance", "y"="Depth.m"))

# plot the smoothed map with value of first axis of the pca
ggplot(di2, aes(x=distance, y=Depth.m)) + geom_raster(aes(fill=value, na.rm=T)) + stat_contour(aes(z=value), colour="white", bins=5, na.rm=TRUE) + scale_fill_gradientn("PCA1", colours=spectral(), na.value=NA) + scale_x_continuous("Distance from shore", expand=c(0,0)) + scale_y_reverse("Depth (m)", expand=c(0,0)) +
opts

#    
# ------------------
# Not too bad, the front is (even more) easy to visualize and appears more vertical than when working on salinity data only.
# Should be tested on other transects.


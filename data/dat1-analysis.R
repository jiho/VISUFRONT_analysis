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
dir <- "/Users/faillettaz/Dropbox/robin/visufront-data/"


library("plyr")
library("ggplot2")
library("stringr")
library("reshape2")
library("grid")
library("gridExtra")
library("akima")
library("fields")

source("lib_zooprocess.R")
source("lib_plot.R")


files <- list.files(dir, full = T)
pid <- read.pid(files[2])



pids <- NA
pids <- adply(files, 1, function(files)
{
pid <- read.pid(files)
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

head(pids)

# set depth BIN 
pids$DepthBin <- round_any(pids$Depth, 1)

# get raw abundance
bio <- ddply(pids, ~Valid+Label+DepthBin, function(x) {sum(na.omit(x$Valid==paste(x$Valid)))})

bio <- rename(bio, c("V1" = "Abund"))

# plot everything
ggplot(bio) + 
geom_point(aes(x = Abund, y = -DepthBin, colour = Label)) + 
geom_path(aes(x = Abund, y = -DepthBin, colour = Label)) + 
facet_grid(.~Valid, scales="free_x")



###         Process Physical Data
### ---------------------------------------------------------

# Read physical data from transect 4
phy <- read.csv("../transects/cross_current_4/isiis.csv", header=T, sep=",")
head(phy)
# delete first line (data from previous transect)
phy <- phy[-1, ]

# select only data from upcasts
d <- phy[which(phy$down.up %in% "up"), ]

# check if seems ok
ggplot(d, aes(x=distanceFromVlfr, y=-Depth.m, colour=Salinity.PPT))  + geom_point() # Yes !






# compute interpolation 
# interpolate all variables
dm <- melt(d, id.vars=c("Depth.m", "distanceFromVlfr"), measure.vars=c("Salinity.PPT"))#, "Temp.C", "Fluoro.volts", "Oxygen.ml.l"))


# First interpolation over a large grid
di <- ddply(dm, ~variable, function(x) {
    x <- na.omit(x)
    xi <- interp.dist(x=x$distanceFromVlfr, y=x$Depth.m, z=x$value, duplicate="mean", x.step=500, y.step=2.5, anisotropy=1300)
})

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




### Second interpolation over a finer grid with spline 
###----------------------------------------------------


# delete NAs from the previous interpolation
i2 <- na.exclude(di)

# Select variables
x <- i2$distance
y <- i2$Depth.m
z <- i2$value


# interpolation with ONLY ONE VARIABLE
# set new step
x.step=1
y.step=1

i <- interp(x=x, y=y, z=z, xo=seq(0, max(x), by=x.step), yo=seq(0, max(y), by=y.step), duplicate="mean", linear=T)


# Second interpolation for all variables
di2 <- ddply(di, ~variable, function(x) {
    x <- na.omit(x)
    xi <- interp.smooth(x=x$distance, y=x$Depth.m, z=x$value, x.step = 0.05, y.step = 0.05)
})


di2 <- rename(di2, c("x"="distance", "y"="Depth.m"))

plots <- dlply(di2, ~variable, function(x) {
        ggplot(x, aes(x=distance, y=-Depth.m)) +
        #geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
        geom_tile(aes(fill=value), na.rm=T) +
        stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
        scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
        })

geom_tile(aes(fill=value), na.rm=T) +
stat_contour(aes(z=value), colour="white", alpha=0.7, bins=6, size=0.2, na.rm=TRUE) +
scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
scale_x_continuous(expand=c(0,0)) +
scale_y_continuous(expand=c(0,0))





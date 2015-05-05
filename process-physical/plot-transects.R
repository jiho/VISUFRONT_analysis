#
#      Read TS record and split it per transect
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#------------------------------------------------------------

# IF NOT IN THE .RPROFILE, SET USER FIRST 
#user <- "faillettaz"

message("Read and process ISIIS hydrological record")

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")

source("lib_process.R")

# dropbox location. change for every user
dropboxloc <- str_c("/Users/", user, "/Dropbox/visufront-data/")

# get all files
tsFiles <- list.files(paste(dropboxloc, "TS", sep=""), pattern="*.tethys", full=TRUE)
ts <- adply(tsFiles, 1, function(file) {
  read.ts(file)
}, .progress="text")
ts <- ts[,-1]


##{ Check the data --------------------------------------------------------

# read coastline to plot the trajectories and check
coast <- read.csv("map/gshhg_coteazur_f.csv")
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast)

# check trajectory
ggplot() + gcoast + geom_point(aes(x=lon, y=lat, colour=-depth), size=1, alpha=0.1, na.rm=T, data = ts) + coord_map()

# remove outliers
ts <- ts[-which(ts$salinity < 31),]

# check T-S diagram
ggplot(ts) + geom_point(aes(x=temperature, y=salinity), size=1.5, alpha=0.1, na.rm=T) + coord_map()

# check time series of all variables
tsm <- melt(ts, id.vars=c("dateTime", "lon", "lat"))
ggplot(tsm) + geom_point(aes(x=dateTime, y=value), size=1.5, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# }


##{ Cut by transect -------------------------------------------------------

# read transects limits
transects <- read.csv("transects.csv", na.strings=c("", "NA"), colClasses=c("character", "POSIXct", "POSIXct"), sep=";")

d <- ddply(transects, ~name, function(x, data) {
	#message(x$name)

  # extract the appropriate portion of the data
  cData <- data[which(data$dateTime > x$dateTimeStart-5 & data$dateTime < x$dateTimeEnd+5),]

  if (nrow(cData) >= 1) {
    # compute distance from first point and from a reference point
    cData$distanceFromStart <- dist.from.start(lat=cData$lat, lon=cData$lon)
    cData$distanceFromVlfr <- dist.from.villefranche(lat=cData$lat, lon=cData$lon)
    cData$distanceFromShore <- dist.from.shore(lat=cData$lat, lon=cData$lon)


    data.frame(cData, transect = x$name)
  }
}, data=ts)


# plot all transects with ID
ggplot() + gcoast + geom_point(aes(x=lon, y=lat, colour=transect), size=3, na.rm=T, data = d) + coord_map()

# plot only selected transects
ggplot() + gcoast + geom_point(aes(x=lon, y=lat, colour=transect), size=3, na.rm=T, data = d[which(substring(d$transect, 1, 5) %in% "cross"), ]) + coord_map()

# plot only selected transects
ggplot() + gcoast + geom_point(aes(x=lon, y=lat, colour=transect), size=1, na.rm=T, data = d[which(d$transect %in% c("cross_current_4", "cross_current_5")), ]) + coord_map()

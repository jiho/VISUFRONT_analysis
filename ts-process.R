#
#      Read TS record and split it per transect
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#------------------------------------------------------------

message("Read and process TS record")

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")

source("lib_process.R")		# function read.ts

# get all files
files <- list.files("TS", pattern="*.tethys", full=TRUE)

ts <- adply(files, 1, function(file) {
  read.ts(file)
}, .progress="text")
ts <- ts[,-1]


##{ Check the data --------------------------------------------------------

# read coastline to plot the trajectories and check
coast <- read.csv("map/gshhg_coteazur_f.csv")
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast)

# check trajectory
ggplot(ts) + gcoast + geom_point(aes(x=lon, y=lat), size=1, alpha=0.1, na.rm=T) + coord_map()
ggplot(ts) + gcoast + geom_point(aes(x=lon, y=lat, colour=-depth), size=1, alpha=0.1, na.rm=T) + coord_map()

# check TS
ts <- ts[-which(ts$salinity < 31),]
ggplot(ts) + geom_point(aes(x=temperature, y=salinity), size=1.5, alpha=0.1, na.rm=T) + coord_map()

# check time series of variables
tsm <- melt(ts, id.vars=c("dateTime", "lon", "lat"))
ggplot(tsm) + geom_point(aes(x=dateTime, y=value), size=1.5, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")
ggplot(tsm) + geom_point(aes(x=lon, y=value), size=1.5, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# }


##{ Cut by transect -------------------------------------------------------

write.csv(ts, file="ts.csv", row.names=FALSE)

transects <- read.csv("transects.csv", na.strings=c("", "NA"), colClasses=c("character", "POSIXct", "POSIXct"))

pdf("ts-transects.pdf")
d_ply(transects, ~name, function(x, data) {
	message(x$name)
	
  # extract the portion of the GPS track
  cData <- data[data$dateTime > x$dateTimeStart-5 & data$dateTime < x$dateTimeEnd+5,]

  if (nrow(cData) >= 1) {  
    # compute distance from first point and from a reference point
    cData$distanceFromStart <- geodDist(lat1=cData$lat, lon1=cData$lon, lat2=cData$lat[1], lon2=cData$lon[1]) * 1.852
    cData$distanceFromVlfr <- geodDist(cData$lat, cData$lon, lat2=43.70528, lon2=7.3118057) * 1.852
    cData$distanceFromShore <- dist.from.shore(cData$lat, cData$lon) / 1.852
  
	  # plot to check
    print(ggplot(cData) + gcoast + geom_point(aes(x=lon, y=lat, colour=distanceFromShore), size=1, alpha=0.1, na.rm=T) + coord_map() + ggtitle(x$name))
  
    # store data file
    dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
    dataName <- deparse(substitute(data))
    write.csv(cData, file=str_c("transects/", x$name, "/", dataName, ".csv"), row.names=FALSE)    
  }
}, data=ts)
dev.off()

# system("open ts-transects.pdf")

# }



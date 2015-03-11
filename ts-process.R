#
#      Read TS record and split it per transect
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#------------------------------------------------------------

source("lib_process.R")
data_dir <- data_dir_path()

library("ggplot2")
library("reshape2")
library("lubridate"))
library("stringr")
library("plyr")
library("dplyr")


##{ Read data --------------------------------------------------------------

# get all files
tsFiles <- list.files(paste(data_dir, "/_raw_/TS", sep=""), pattern="*.tethys", full=TRUE)
ts <- ldply(tsFiles, function(file) {
  read.ts(file)
}, .progress="text")

# }

##{ Check the data --------------------------------------------------------

# read coastline to plot the trajectories and check
coast <- read.csv("map/gshhg_coteazur_f.csv")
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast)

# check trajectory
ggplot(ts) + gcoast + geom_point(aes(x=lon, y=lat, colour=-depth), size=1, alpha=0.1, na.rm=T) + coord_map()

# remove outliers
ts <- ts[-which(ts$salinity < 31),]

# check T-S diagram
ggplot(ts) + geom_point(aes(x=temperature, y=salinity), size=1.5, alpha=0.1, na.rm=T)

# check time series of all variables
tsm <- melt(ts, id.vars=c("dateTime", "lon", "lat"))
ggplot(tsm) + geom_point(aes(x=dateTime, y=value), size=1.5, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# }


##{ Cut by transect -------------------------------------------------------

# read transects limits
transects <- read.csv(str_c(data_dir, "/transects.csv"), na.strings=c("", "NA"))
transects$dateTimeStart <- ymd_hms(transects$dateTimeStart, tz="Europe/Paris")
transects$dateTimeEnd <- ymd_hms(transects$dateTimeEnd, tz="Europe/Paris")

# add transect name to the file


pdf("ts-transects.pdf")
ts_in_transects <- ddply(transects, ~name, function(x, data) {
	message(x$name)

  # extract the appropriate portion of the data
  cData <- filter(data, dateTime > x$dateTimeStart-5, data$dateTime < x$dateTimeEnd+5)

  if (nrow(cData) >= 1) {
    # compute distance from first point and from a reference point
    cData$distance_from_start <- dist.from.start(lat=cData$lat, lon=cData$lon)
    cData$distance_from_vlfr <- dist.from.villefranche(lat=cData$lat, lon=cData$lon)
    cData$distance_from_shore <- dist.from.shore(lat=cData$lat, lon=cData$lon)

	  # plot to check
    print(ggplot(cData) + gcoast + geom_point(aes(x=lon, y=lat, colour=distance_from_shore), size=1, alpha=0.3, na.rm=T) + coord_map() + ggtitle(x$name))

    # store data file
    output_dir <- str_c(data_dir, "/transects/", x$name)
    dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
    write.csv(cData, file=str_c(output_dir, "/ts.csv"), row.names=FALSE)
  }

	return(cData)
}, data=ts)
dev.off()

# system("open ts-transects.pdf")

# write the complete record
write.csv(ts, file=str_c(data_dir, "/ts_all.csv"), row.names=FALSE)
write.csv(ts_in_transects, file=str_c(data_dir, "/ts_in_transects.csv"), row.names=FALSE)

# }



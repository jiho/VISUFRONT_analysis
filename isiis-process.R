#
#      Read ISIIS physical data and split it per transect
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#------------------------------------------------------------

message("Read and process ISIIS hydrological record")

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")

source("lib_process.R")

# dropbox location. change for every user
dropboxloc <- "/Users/jessicaluo/Dropbox/"

# get data
hydroFiles <- list.files(paste(dropboxloc, "visufront-data/ISIIShydro", sep=""), pattern=glob2rx("ISIIS*.txt"), full=T)
d <- adply(hydroFiles, 1, function(file) {
	read.isiis(file)
}, .progress="text")
d <- d[,-1]


##{ Cleanup data ----------------------------------------------------------

# remove erroneous values
d$Temp.C[d$Temp.C <= 0] <- NA
d$Salinity.PPT[d$Salinity.PPT <= 10] <- NA
d$Fluoro.volts[d$Fluoro.volts <= 0.01] <- NA
d$Oxygen.ml.l[d$Oxygen.ml.l <= 0] <- NA
d$Oxygen.ml.l[d$Oxygen.ml.l >= 8] <- NA

# check TS diagram
ggplot(d) + geom_point(aes(x=Temp.C, y=Salinity.PPT), size=1, alpha=0.1, na.rm=T)

# check all variables
dm <- melt(d, id.vars=c("dateTimeMsec", "Pressure.dbar", "Depth.m"))
ggplot(dm) + geom_point(aes(x=dateTimeMsec, y=value), size=1, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# }


##{ Add location data -----------------------------------------------------

# read TS record, which contains GPS location
ts <- read.csv("ts.csv", stringsAsFactors=FALSE)
ts$dateTime <- as.POSIXct(ts$dateTime, tz="GMT") # we are not in GMT but this is added to get rid of an error with the join later
# round ISIIS time to the second, to match with the ship's GPS
d$dateTime <- round(d$dateTimeMsec)
# get lat-lon from the TS record

# TODO: error here "Error in as.POSIXct.POSIXlt(what, tz = tzone) : invalid 'tz' value"
# Robin? Any ideas to fix this?
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")
d <- join(d, ts[,c("dateTime", "lat", "lon")], by="dateTime")

# interpolate GPS data
sum(is.na(d$lat))
sum(is.na(d$lon))
d$lat <- approx(x=as.numeric(d$dateTimeMsec), y=d$lat, xo=as.numeric(d$dateTimeMsec))$y # another error here "Error in xy.coords(x, y) : 'x' and 'y' lengths differ"
d$lon <- approx(x=as.numeric(d$dateTimeMsec), y=d$lon, xo=as.numeric(d$dateTimeMsec))$y # error "Error in xy.coords(x, y) : 'x' and 'y' lengths differ"
sum(is.na(d$lat))
sum(is.na(d$lon))

# use salinity, temperature and pressure to compute seawater density using UNESCO formulation
d$Density <- swRho(d$Salinity.PPT, d$Temp.C, d$Pressure.dbar, eos="unesco")

# }


##{ Cut by transect -------------------------------------------------------

# write the full record
isiis <- d
write.csv(isiis, file="isiis.csv", row.names=FALSE)

# read transects limits
transects <- read.csv("transects.csv", na.strings=c("", "NA"), colClasses=c("character", "POSIXct", "POSIXct"))

pdf("isiis-transects.pdf")
d_ply(transects, ~name, function(x, data) {
	message(x$name)

  # extract the appropriate portion of the data
  cData <- data[which(data$dateTime > x$dateTimeStart-5 & data$dateTime < x$dateTimeEnd+5),]

  if (nrow(cData) >= 1) {
    # compute distance from first point and from a reference point
    cData$distanceFromStart <- dist.from.start(lat=cData$lat, lon=cData$lon)
    cData$distanceFromVlfr <- dist.from.villefranche(lat=cData$lat, lon=cData$lon)
    cData$distanceFromShore <- dist.from.shore(lat=cData$lat, lon=cData$lon)

    # detect up and down casts
		casts <- detect.casts(cData$Depth.m)
		cData <- cbind(cData, casts)

    # plot to check
    print(ggplot(cData) + geom_point(aes(x=distanceFromShore, y=-Depth.m, colour=down.up), na.rm=T) + ggtitle(x$name))

    # store data file
    dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
    dataName <- deparse(substitute(data))
    write.csv(cData, file=str_c("transects/", x$name, "/", dataName, ".csv"), row.names=FALSE)
  }

	return(invisible(NULL))
}, data=isiis)
dev.off()

# }

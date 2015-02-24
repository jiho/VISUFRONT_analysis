#
#      Read ISIIS physical data and split it per transect
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

data <- "~/Dropbox/visufront-data/"

library("lubridate")
library("stringr")
library("ggplot2")
library("reshape2")
library("plyr")
library("dplyr")

source("lib_process.R")


##{ Read ISIIS hydro data --------------------------------------------------

# get data
hydroFiles <- list.files(str_c(data, "ISIIShydro"), pattern=glob2rx("ISIIS*.txt"), full=TRUE)
d <- ldply(hydroFiles, function(file) {
	read.isiis(file)
}, .progress="text")

# }


##{ Cleanup data ----------------------------------------------------------

# remove erroneous values
d$Temp.C[d$Temp.C <= 0] <- NA
d$Salinity.PPT[d$Salinity.PPT <= 10] <- NA
d$Fluoro.volts[d$Fluoro.volts <= 0.01] <- NA
d$Oxygen.ml.l[d$Oxygen.ml.l <= 0] <- NA
d$Oxygen.ml.l[d$Oxygen.ml.l >= 8] <- NA

# use salinity, temperature and pressure to compute seawater density using UNESCO formula
d$Density <- swRho(d$Salinity.PPT, d$Temp.C, d$Pressure.dbar, eos="unesco")

# # check TS diagram
# ggplot(d) + geom_point(aes(x=Temp.C, y=Salinity.PPT), size=1, alpha=0.1, na.rm=T)
#
# # check all variables
# dm <- melt(d, id.vars=c("dateTimeMsec", "Pressure.dbar", "Depth.m"))
# ggplot(dm) + geom_point(aes(x=dateTimeMsec, y=value), size=1, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# }


##{ Add lat/lon -----------------------------------------------------------

# read TS record, which contains GPS location
ts <- read.csv(str_c(data, "TS/ts.csv"), stringsAsFactors=FALSE)
ts$dateTime <- ymd_hms(ts$dateTime)

# round ISIIS time to the second, to match with the ship's GPS
d$dateTime <- round_any(d$dateTimeMsec, 1)
# NB: round turns this into a POSIXlt object

# get lat-lon from the TS record
d <- left_join(d, select(ts, dateTime, lat, lon), by="dateTime")

# interpolate GPS data
sum(is.na(d$lat))
sum(is.na(d$lon))
d$lat <- approx(x=as.numeric(d$dateTimeMsec), y=d$lat, xo=as.numeric(d$dateTimeMsec))$y
d$lon <- approx(x=as.numeric(d$dateTimeMsec), y=d$lon, xo=as.numeric(d$dateTimeMsec))$y 
sum(is.na(d$lat))
sum(is.na(d$lon))

# }


##{ Cut by transect -------------------------------------------------------

# read transects limits
transects <- read.csv(str_c(data, "transects.csv"), na.strings=c("", "NA"))
transects$dateTimeStart <- ymd_hms(transects$dateTimeStart)
transects$dateTimeEnd <- ymd_hms(transects$dateTimeEnd)

pdf("isiis-transects.pdf")
isiis_in_transects <- ddply(transects, ~name, function(x, data) {
	message(x$name)

  # extract the appropriate portion of the data
  cData <- filter(data, dateTime > x$dateTimeStart-5, data$dateTime < x$dateTimeEnd+5)

  if (nrow(cData) >= 1) {
    # compute distance from first point and from a reference point
    cData$dist_from_start <- dist.from.start(lat=cData$lat, lon=cData$lon)
    cData$dist_from_vlfr <- dist.from.villefranche(lat=cData$lat, lon=cData$lon)
    cData$dist_from_shore <- dist.from.shore(lat=cData$lat, lon=cData$lon)

    # detect up and down casts
    casts <- detect.casts(cData$Depth.m)
    cData <- cbind(cData, casts)

    # plot to check
    print(ggplot(cData) + geom_point(aes(x=dist_from_shore, y=-Depth.m, colour=down.up), na.rm=T) + ggtitle(x$name))

    # remove incorrect data: above 30m in downcasts
    cData <- filter(cData, !(cData$down.up %in% "down" & cData$Depth.m <= 30))
    # ggplot(e, aes(x=distanceFromVlfr, y=-Depth.m, colour=down.up)) + geom_point()
    # ggplot(eC, aes(x=distanceFromVlfr, y=-Depth.m, colour=Salinity.PPT))  + geom_point() + scale_colour_spectral()

    # store data file
    dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
    write.csv(cData, file=str_c("transects/", x$name, "/isiis.csv"), row.names=FALSE)
  }

	return(cData)
}, data=d)
dev.off()

# write the full record
write.csv(d, file="isiis.csv", row.names=FALSE)
write.csv(isiis_in_transects, file="isiis_in_transects.csv", row.names=FALSE)

# }

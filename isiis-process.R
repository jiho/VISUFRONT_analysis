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
library("oce")

source("lib_process.R")


##{ Read ISIIS hydro data --------------------------------------------------

# get data
hydroFiles <- list.files(str_c(data, "ISIIShydro"), pattern=glob2rx("ISIIS*.txt"), full=TRUE)
d <- ldply(hydroFiles, function(file) {
	read.isiis(file)
}, .progress="text")

options(digits.secs=2)
d$dateTimeMsec <- ymd_hms(d$dateTimeMsec)

# }


##{ Cleanup data ----------------------------------------------------------

# remove an small back and forth trajectory in cross-front transect 5
# ISIIS camera failure
d <- filter(d, !(dateTimeMsec > ymd_hms("2013-07-25 08:11:30") & dateTimeMsec < ymd_hms("2013-07-25 10:32:11")))

# remove erroneous values
d$Temp.C[d$Temp.C <= 0] <- NA
d$Salinity.PPT[d$Salinity.PPT <= 10] <- NA
d$Fluoro.volts[d$Fluoro.volts <= 0.01] <- NA
d$Oxygen.ml.l[d$Oxygen.ml.l <= 0] <- NA
d$Oxygen.ml.l[d$Oxygen.ml.l >= 8] <- NA

# use salinity, temperature and pressure to compute seawater density using UNESCO formula
d$Density <- swRho(d$Salinity.PPT, d$Temp.C, d$Pressure.dbar, eos="unesco")

# identify the role of variables
vars <- c("Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "Density")
id_vars <- c("dateTimeMsec", "Pressure.dbar", "Depth.m")

# # check TS diagram
# ggplot(d) + geom_point(aes(x=Temp.C, y=Salinity.PPT), size=1, alpha=0.1, na.rm=T)
#
# # check all variables
# dm <- melt(d, id.vars=id_vars)
# ggplot(dm) + geom_point(aes(x=dateTimeMsec, y=value), size=1, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# }


##{ Add lat/lon -----------------------------------------------------------

# read TS record, which contains GPS location
ts <- read.csv(str_c(data, "TS/ts.csv"), stringsAsFactors=FALSE)
ts$dateTime <- ymd_hms(ts$dateTime)

# round ISIIS time to the second, to match with the ship's GPS
d$dateTime <- round_any(d$dateTimeMsec, 1)
# NB: round turns this into a POSIXlt object, hence the use of round_any

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
isiis_in_transects <- ddply(transects, ~name, function(x, d) {
	message(x$name)

  # extract the appropriate portion of the data
  df <- filter(d, dateTime > x$dateTimeStart-5, d$dateTime < x$dateTimeEnd+5)

  if (nrow(df) >= 1) {
    # compute distance from first point and from a reference point
    df$dist_from_start <- dist.from.start(lat=df$lat, lon=df$lon)
    df$dist_from_vlfr <- dist.from.villefranche(lat=df$lat, lon=df$lon)
    df$dist_from_shore <- dist.from.shore(lat=df$lat, lon=df$lon)

    # detect up and down casts
    casts <- detect.casts(df$Depth.m)
    df <- cbind(df, casts)

    # plot to check
    print(ggplot(df) + geom_point(aes(x=dist_from_shore, y=-Depth.m, colour=down.up), na.rm=T) + ggtitle(x$name))

    # remove incorrect data
    # salinity above 30m in downcasts
    df[which(df$down.up %in% "down" & df$Depth.m <= 30), c("Salinity.PPT", "Density")] <- NA
    # ggplot(df, aes(x=dist_from_start, y=-Depth.m, colour=Salinity.PPT))  + geom_point() + scale_color_spectral()
    # and irradiance
    df <- select(df, -Irrandiance.UE.cm)

    # store data file
    dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
    write.csv(df, file=str_c("transects/", x$name, "/isiis.csv"), row.names=FALSE)
  }

	return(df)
}, d=d)
dev.off()

# write the full record
write.csv(d, file="isiis.csv", row.names=FALSE)
write.csv(isiis_in_transects, file="isiis_in_transects.csv", row.names=FALSE)

# }

#
#      Plot drifters trajectories
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("gdata")
library("stringr")
library("plyr")
library("ggplot2")
library("oce")

source("lib_process.R")

#--------------------------------------------------------------------------

outliers <- function(x, method=c("hampel","g","bonferroni","custom"), factor=5.2) {
  #
  #	Return the indices of outliers in x, according to:
  #		. Davies and Gather, The identification of multiple outliers, JASA 88 (1993), 782-801. for methods hampel, g and custom
  #		. outlier.test in package car for method bonferroni
  #	In method custom, the higher the factor the less sensible the detection of outliers
  #
	method = match.arg(method)

	if (method=="bonferroni") {
		suppressPackageStartupMessages(require("car"))
		return(as.numeric(outlierTest(lm(x~1))$obs))
	} else {
		if (method=="hampel") {
			factor = 5.2
		} else if (method=="g") {
    		n = length(x)
			if (n%%2==0) {
				factor = 2.906+11.99*(n-6)^-0.5651
			} else {
				factor = 2.906+12.99*(n-5)^-0.5781
			}
		} else if (method=="custom") {
			factor = factor
		}
		return(which(abs(x-median(x, na.rm=TRUE))>(factor*mad(x, na.rm=TRUE))))
	}
}

#--------------------------------------------------------------------------

# read coastline
coast <- read.csv("map/gshhg_coteazur_f.csv")
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast)

# read drifter trajectories
d <- read.xls("drifters.xls", na.strings="", stringsAsFactors=FALSE)
b <- read.xls("boa.xls", na.strings="", stringsAsFactors=FALSE)

# read ship trajectory from ts
s1 <- read.ts(str_c("/Volumes/donnees/", format(Sys.Date()-1, "%Y%m%d"), ".tethys"))
s2 <- read.ts(str_c("/Volumes/donnees/", format(Sys.Date(), "%Y%m%d"), ".tethys"))
s <- rbind(s1,s2)

# fill date
for (i in 1:nrow(d)) {
  if (is.na(d$date[i])) {
    d$date[i] <- d$date[i-1]
  }
}
for (i in 1:nrow(b)) {
  if (is.na(b$date[i])) {
    b$date[i] <- b$date[i-1]
  }
}

# compute lat long for boa
latBits <- str_split_fixed(b$lat, fixed("."), 2)
b$lat <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
b$lat <- 43 + b$lat/60
lonBits <- str_split_fixed(b$lon, fixed("."), 2)
b$lon <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
b$lon <- 7 + b$lon/60

# combine drifters and boa
d <- rbind(d,b)

# sort by date and time
d$dateTime <- as.POSIXct(str_c(d$date, d$time, sep=" ")) + 2 * 3600
d <- arrange(d, dateTime, unit)

# recode the units to fit the ship's notation
d$sunit <- as.character(d$unit)
d$sunit <- factor(d$sunit, levels=c("20","30","provbio","boa"), labels=c("C1", "C2", "provbio", "C3"))
d$unit <- str_c(d$sunit, d$unit, sep=" - ")

# add a record number per float
d <- ddply(d, ~unit, function(x) {
  x$nb <- 1:nrow(x)
  return(x)
})

# add a time label
d$timeLabel <- format(d$dateTime, format="%H:%M")

# ggplot(d, aes(x=dateTime, y=lat, colour=unit)) + geom_path(na.rm=T)
# ggplot(d, aes(x=dateTime, y=lon, colour=unit)) + geom_path(na.rm=T)

# remove outliers in lat and lon
d$lon[d$lon > 8 | d$lon < 7] <- NA
d$lon[outliers(d$lon)] <- NA
d$lat[outliers(d$lat)] <- NA

# remove faulty drifters
# d <- d[!d$unit %in% c("20", "30"),]

# interpolate drifter position to the current position of the ship
di <- ddply(d, ~unit, function(x, ship) {
	now <- ship$dateTime[nrow(ship)]
	iLat <- spline(x=x$dateTime, y=x$lat, xout=now, method="natural")$y
	iLon <- spline(x=x$dateTime, y=x$lon, xout=now, method="natural")$y
	return(data.frame(dateTime=now, lat=iLat, lon=iLon))
}, ship=s)

# compute lat-lon compatible with the ship's system
d$latMin <- (d$lat - floor(d$lat)) * 60
d$lonMin <- (d$lon - floor(d$lon)) * 60
di$latMin <- (di$lat - floor(di$lat)) * 60
di$lonMin <- (di$lon - floor(di$lon)) * 60


# # get ship data corresponding to the drifters trajectories
# startTime <- min(d$dateTime)
# stopTime <- max(d$dateTime)
# s <- s[s$dateTime >= startTime & s$dateTime <= stopTime,]

# set zoom limits
limits <- ddply(d[d$unit != "provbio - provbio",], ~unit, function(x) {
  x <- tail(x)
  data.frame(min(x$lat), max(x$lat), min(x$lon), max(x$lon))
})
limits <- with(limits, c(min(min.x.lat.), max(max.x.lat.), min(min.x.lon.), max(max.x.lon.)))
pad <- 0.02
limits <- limits + c(-pad, pad, -pad, pad)

pdf("drifters.pdf", width=10, height=7)
print(
ggplot(d, aes(x=lon, y=lat)) +
  coord_map() +
  gcoast +
  # ship track
  geom_path(size=0.1, na.rm=T, data=s, alpha=0.8) +
  # drifters track
  geom_path(aes(colour=unit), size=0.1, na.rm=T, alpha=1) +
  geom_point(aes(colour=unit), size=0.3, na.rm=T, alpha=1) +
  geom_text(aes(colour=unit, label=timeLabel), size=0.5, hjust=1.1, vjust=-0.3, na.rm=T) +
	# interpolated current position
  geom_point(aes(colour=unit), size=0.5, na.rm=T, alpha=1, data=di) +
  scale_x_continuous(expand=c(0.01,0.01))
  # xlim(limits[3], limits[4]) + ylim(limits[1], limits[2])
)
dev.off()
# print(ggplot(d, aes(x=lon, y=lat)) + gcoast + geom_path(aes(colour=unit), na.rm=T) + geom_text(aes(colour=unit, label=nb), size=3, hjust=0, vjust=0.5, na.rm=T) + facet_wrap(~unit, scales="free"))

# display the last positions
# message("Last position")
# print(ddply(d, ~unit, function(x) {
#   tail(x, 1)[,c("unit", "dateTime", "latMin", "lonMin")]
# }))

# estimate speed
message("Speed")
speed <- ddply(d, ~unit, function(x) {
  x <- tail(x)
  dist <- geodDist(lat1=x$lat, lon1=x$lon, alongPath=TRUE)
  time <- as.numeric(difftime(max(x$dateTime), min(x$dateTime), units="secs"))
  speedCMS <- dist[length(dist)] * 1000 * 100 / time
  speedKNTS <- dist[length(dist)] / 1.852 / (time / 3600)
  data.frame(speedCMS, speedKNTS)
})
print(speed)

message("Predicted current position")
print(di[,c("unit", "dateTime", "latMin", "lonMin")])







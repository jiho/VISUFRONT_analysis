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

# read drifter trajectories
d <- read.xls("drifters/drifters.xls", na.strings="", stringsAsFactors=FALSE)
b <- read.xls("drifters/boa.xls", na.strings="", stringsAsFactors=FALSE)
f <- read.xls("drifters/float.xls", na.strings="", stringsAsFactors=FALSE)

# read ship trajectory from ts
s1 <- read.ts("TS/20130727.tethys")
s2 <- read.ts("TS/20130728.tethys")
s <- rbind(s1,s2)

# compute lat and lon for boa
latBits <- str_split_fixed(b$lat, fixed("."), 2)
b$lat <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
b$lat <- 43 + b$lat/60
lonBits <- str_split_fixed(b$lon, fixed("."), 2)
b$lon <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
b$lon <- 7 + b$lon/60

# combine all data
d <- rbind(d,b,f)

# sort by date and time
d$dateTime <- as.POSIXct(str_c(d$date, d$time, sep=" ")) + 2 * 3600
d <- arrange(d, dateTime, unit)

# recode the name of the drifters to fit the ship's notation
sunit <- as.character(d$unit)
sunit <- factor(sunit, levels=c("20","30","provbio","boa"), labels=c("C1", "C2", "provbio", "C3"))
d$unit <- str_c(sunit, d$unit, sep=" - ")

# add a record number per float
d <- ddply(d, ~unit, function(x) {
  x$nb <- 1:nrow(x)
  return(x)
})

# add a time label to show on the plot
d$timeLabel <- format(d$dateTime, format="%H:%M")

# interpolate drifter position at the present time or in the future
# 	ship 	is the ship track form the TS
#		lag		is a time lag in hours
di <- ddply(d, ~unit, function(x, ship, lag=0) {
	target <- Sys.time() + (lag * 3600)
	iLat <- spline(x=x$dateTime, y=x$lat, xout=target, method="natural")$y
	iLon <- spline(x=x$dateTime, y=x$lon, xout=target, method="natural")$y
	return(data.frame(dateTime=target, lat=iLat, lon=iLon))
}, ship=s)

# compute lat-lon compatible with the ship's system
d$latMin <- (d$lat - floor(d$lat)) * 60
d$lonMin <- (d$lon - floor(d$lon)) * 60
di$latMin <- (di$lat - floor(di$lat)) * 60
di$lonMin <- (di$lon - floor(di$lon)) * 60

# read and pre-plot coastline
coast <- read.csv("map/gshhg_coteazur_f.csv")
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast)

# # get ship data corresponding to the drifters trajectories
startTime <- min(d$dateTime)
stopTime <- max(d$dateTime[d$unit=="C1 - 20"])
s <- s[s$dateTime >= startTime & s$dateTime <= stopTime,]

pdf("drifters.pdf", width=10, height=7)
print(
ggplot(d, aes(x=lon, y=lat)) +
  coord_map() +
  gcoast +
  # ship track
  geom_path(size=0.5, na.rm=T, data=s, alpha=0.8) +
  # drifters track
  geom_path(aes(colour=unit), size=0.5, na.rm=T, alpha=1) +
  geom_point(aes(colour=unit), size=1, na.rm=T, alpha=1) +
  # geom_text(aes(colour=unit, label=timeLabel), size=0.5, hjust=1.1, vjust=-0.3, na.rm=T) +
	# interpolated position
  geom_point(aes(colour=unit), size=1.5, na.rm=T, alpha=1, data=di)
  # scale_x_continuous(expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0))
)
dev.off()

# Display information messages

# message("Last position")
# print(ddply(d, ~unit, function(x) {
#   tail(x, 1)[,c("unit", "dateTime", "latMin", "lonMin")]
# }))

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

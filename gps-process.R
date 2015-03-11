#
#      Read and cut GPS data
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

source("lib_process.R")
data_dir <- data_dir_path()

library("plyr")
library("stringr")
library("oce")
library("dplyr")
library("ggplot2")
library("lubridate")
library("discuss")
library("circular")


# Concatenate all GPS files
gpsFiles <- list.files(str_c(data_dir, "/_raw_/gps"), pattern=glob2rx("gps*"), full=TRUE)
gps <- alply(gpsFiles, 1, function(file) {
  # read the content of the file
  content <- scan(file, what="character", sep="\n", fileEncoding="ISO-8859-1", encoding="UTF-8")
  return(content)
})
gps <- unlist(gps)
head(gps)

# only keep GPZDA (date) and GPGLL (lat, lon) lines
gpsD <- gps[str_detect(gps, fixed("$GPGGA"))]
gpsL <- gps[str_detect(gps, fixed("$GPGLL"))]


# # fix problems
# gpsD <- str_replace(gpsD, fixed("07<2013"), "07,2013")
# gpsD <- str_replace(gpsD, fixed("+0\u0082àR¢j"), "+00,00")
# gpsD <- str_replace(gpsD,fixed("$GPGGA,072104.00,4325.99442,N,00747.22920,E,1,09,0.9,17.39,M"), "")
gpsD <- gpsD[-85986]
gpsD <- gpsD[-162243]
gpsD <- gpsD[-273354]
gpsD <- gpsD[-275104]
gpsD <- gpsD[-326743]
gpsD <- gpsD[-523406]

gpsL <- gpsL[-298594]
gpsL <- gpsL[-377905]

# Read concatenated data
read.gps <- function(x) {
  temp <- tempfile()
  cat(x, file=temp, sep="\n")
  x <- read.table(temp, sep=",", header=FALSE, stringsAsFactors=FALSE, colClasses="character")
  file.remove(temp)
  return(x)
}

gpsD <- read.gps(gpsD)
gpsL <- read.gps(gpsL)
head(gpsD)
head(gpsL)

# Extract coordinates of gpsL
# the current format is DDMM.MMM
deglat <- as.numeric(str_sub(gpsL$V2, 1, 2))
minlat <- as.numeric(str_sub(gpsL$V2, 3))
deglon <- as.numeric(str_sub(gpsL$V4, 1, 3))
unique(deglon)
deglon[which(is.na(deglon))] <- 7 # Checked manually and all lon should be 7°

minlon <- as.numeric(str_sub(gpsL$V4, 4))
lat <- deglat + minlat/60
lon <- deglon + minlon/60

which(is.na(lon)) == which(is.na(lat)) 

# Extract time
day <- str_sub(gpsL$V1, 10, 11)
month <- "07"
year <- "2013"

hour <- str_sub(gpsL$V6, 1, 2)
min <- str_sub(gpsL$V6, 3, 4)
sec <- str_sub(gpsL$V6, 5, 6)

dateTime <- str_c(str_c(year, month, day, sep="-"), " ", str_c(hour, min, sec, sep=":"))
dateTimeUTC <- ymd_hms(dateTime)


# compensate for time difference with computer time
# dateTime <- dateTimeUTC + 2 * 3600 - 7


# Cleanup unused columns
gpsLf <- data.frame(lat, lon, dateTimeUTC)
gpsLf <- filter(gpsLf, !is.na(lon), !is.na(lat))




# Do the same for gpsD because some data are missing in gpsL
# the current format is DDMM.MMM
deglat <- as.numeric(str_sub(gpsD$V3, 1, 2))
# deglat[which(is.na(deglat))] <- 43
minlat <- as.numeric(str_sub(gpsD$V3, 3))
deglon <- as.numeric(str_sub(gpsD$V5, 1, 3))
# deglon[which(is.na(deglon))] <- 7
minlon <- as.numeric(str_sub(gpsD$V5, 4))
lat <- deglat + minlat/60
lon <- deglon + minlon/60

# check for NAs
which(is.na(lat)) == which(is.na(lon))
length(lat)
length(lon)

# Extract time
day <- str_sub(gpsD$V1, 10, 11)
month <- "07"
year <- "2013"

hour <- str_sub(gpsD$V1, 13, 14)
min <- str_sub(gpsD$V1, 16, 17)
sec <- str_sub(gpsD$V1, 19, 20)

dateTime <- str_c(str_c(year, month, day, sep="-"), " ", str_c(hour, min, sec, sep=":"))
dateTimeUTC <- ymd_hms(dateTime) - 2 * 3600# - 7

# compensate for time difference with computer time
# dateTime <- dateTimeUTC + 2 * 3600 - 7


# Cleanup unused columns
gpsDf <- data.frame(lat, lon, dateTimeUTC)
head(gpsDf)
tail(gpsDf)

gpsDf <- filter(gpsDf, !is.na(lon))


# Cut into transects
# transects <- read.csv("transects.csv", na.strings="", sep = ";", colClasses=c("character", "POSIXct", "POSIXct"))
#
# d_ply(transects, ~name, function(x, gps) {
#   # extract the portion of the GPS track
#   cGps <- gps[gps$dateTime > x$dateTimeStart-5 & gps$dateTime < x$dateTimeEnd+5,]
#
#   if (nrow(cGps) >= 1) {
#     # compute distance from first point and from a reference point
#     cGps$distanceFromStart <- geodDist(lat1=cGps$lat, lon1=cGps$lon, lat2=cGps$lat[1], lon2=cGps$lon[1]) * 1.852
#     cGps$distanceFromVlfr <- geodDist(cGps$lat, cGps$lon, lat2=43.70528, lon2=7.3118057) * 1.852
#
#     # store data file
#     dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
#     write.csv(cGps, file=str_c("transects/", x$name, "/gps.csv"), row.names=FALSE)
#   }
#
# }, gps=gps)





# Compute the bearing of the boat to retrieve the ADCP data
# --------------------------------------------------------------------

# Subset data to work on fewer first
gpsLag <- filter(gpsDf, dateTimeUTC > ymd_hms("2013-07-26 21:30:00") & dateTimeUTC < ymd_hms("2013-07-28 05:25:00"))
range(gpsLag$dateTimeUTC)

# plot it
ggplot() + geom_path(data= gpsLag, aes(x = lon, y = lat, colour = dateTimeUTC))
# SOME DATA ARE MISSING AFTER THE FIRST LAGRANGIAN TRANSECT. We used to have these data

# Find which data are missing
# filter(gpsLag, dateTimeUTC > as.POSIXct("2013-07-27 03:40:00") & dateTimeUTC < as.POSIXct("2013-07-27 09:10:00"))
# NO DATA BETWEEN 3:42 AM to 9:05 AM on the 27th
# 43.55036 7.586471 2013-07-27 03:42:55
# 43.60356 7.417247 2013-07-27 09:05:40


# Check very fine scale path 
p1 <- ggplot(data = filter(filter(gpsLf, dateTimeUTC > ymd_hms("2013-07-26 21:30:00") & dateTimeUTC < ymd_hms("2013-07-28 05:25:00")), dateTimeUTC > ymd_hms("2013-07-26 22:00:00") & dateTimeUTC < ymd_hms("2013-07-26 22:01:00")), aes(x = lon, y = lat)) + geom_path() + geom_point() + ggtitle("$GPGLL")
p2 <- ggplot(data = filter(filter(gpsDf, dateTimeUTC > ymd_hms("2013-07-26 21:30:00") & dateTimeUTC < ymd_hms("2013-07-28 05:25:00")), dateTimeUTC > ymd_hms("2013-07-26 22:00:00") & dateTimeUTC < ymd_hms("2013-07-26 22:01:00")), aes(x = lon, y = lat)) + geom_path() + geom_point() + ggtitle("$GPGGA")


# pdf(file = "plot/gps_1m.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
# dev.off()
# Very straight, this is good. 




# Read TS file in which we have the GPS data every 15s

# Read thermosalinometer from Tethys
# get all files
tsFiles <- list.files("/Users/faillettaz/Dropbox/visufront-data/TS", pattern="*.tethys", full=TRUE)

ts <- adply(tsFiles, 1, function(file) {
	
	# file <- tsFiles[1]
  # get column names from "notice_daufin-TS.pdf" page 12
  colNames <- c(
      "errorCode", "timeUTC", "lat", "lon", "PDOP", "nbSatellites", "sep",
      "timeUTC_human", "latDegMin", "latNS", "lonDegMin", "lonEW", "gpsQuality", "nbSatellites", "HDOP", "sep",
      "atmPressure", "tempAir", "humidity", "windDirection", "windSpeed", "sep",
      "bearing", "shipSpeed", "sep",
      "bearingTrue", "shipSpeedTrue", "sep",
      "depth", "sep",
      "tempCell", "conductivity", "temperature", "salinity", "fluorometry",
      "empty", "empty", "empty", "empty"
  )
  # read data
  d <- read.table(file, sep="\t", header=FALSE, col.names=colNames)

  # create date+time in R format
  year <- str_sub(file, -15, -12)
  month <- str_sub(file, -11, -10)
  day <- str_sub(file, -9, -8)
  midnight <- ymd_hms(str_c(year, "-", month, "-", day, " 00:00:01"))
  d$dateTimeUTC <- midnight + d$timeUTC
  # convert into local time
  # d$dateTime <- d$dateTime + 2*3600

  # fill the gaps in GPS coordinates when necessary
  if (any(is.na(d$lat))) {
    d$lat <- approx(x=d$dateTimeUTC, y=d$lat, xo=d$dateTimeUTC)$y
    d$lon <- approx(x=d$dateTimeUTC, y=d$lon, xo=d$dateTimeUTC)$y
  }
  # select interesting data columns
  return(d)
  
}, .progress = "text")
ts <- ts[,-1]
head(ts)


# # We have access to the bearing of the boat every 15s


# tmp <- sample_frac(ts, 0.05)
# ggplot() + geom_bar(aes(x = bearing), data = tmp, binwidth = 5) + polar()
# # It looks relevent


# Extract data from the lagrangian experiment
tsLag <- filter(ts, dateTimeUTC > ymd_hms("2013-07-26 17:07:00") & dateTimeUTC < ymd_hms("2013-07-28 02:10:00"))
# For some reasons, the time are wrong
ggplot() + geom_path(data= tsLag, aes(x = lon, y = lat))

# Check bearings distribution
ggplot() + geom_bar(aes(x = bearing), data = tsLag, binwidth = 5) + polar()


# Check very fine scale path 
ggplot(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 21:33:00") & dateTimeUTC < ymd_hms("2013-07-26 21:35:00")), aes(x = lon, y = lat)) + geom_path() + geom_point()
# Very straight, this is good.



# Compare gps from TS and gps GPGGA
ggplot() + 
# TS data
geom_path(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:09:50")), aes(x = lon, y = lat), colour = "red") + geom_point(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:09:50")), aes(x = lon, y = lat), colour = "red") +
# GPS GPGGA data
geom_path(data = filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:09:50")), aes(x = lon, y = lat)) + geom_point(data = filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:09:50")), aes(x = lon, y = lat)) + coord_cartesian(xlim = c(7.56, 7.562), ylim = c(43.632, 43.634))
# ggsave(file = "plot/gps_TS_GPGAA.pdf")




# Compute bearing of from the GPS data for the 2 min period
# --------------------------------------------------------------------
gps2m <- filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:09:50"))
gps2m$angle.rad <- NA


# Compute bearing between each point and the next one (1s)
for (i in 1:(nrow(gps2m)-1)) {
    a <- geodDist(gps2m$lat[i+1], gps2m$lon[i+1], gps2m$lat[i], gps2m$lon[i])
    b <- geodDist(gps2m$lat[i], gps2m$lon[i], gps2m$lat[i+1], gps2m$lon[i])
    # check if positive values and correct if 
    if (gps2m$lon[i+1] < gps2m$lon[i]) {
    a <- -a
    }
    if (gps2m$lat[i+1] < gps2m$lat[i]) {
    b <- -b
    }
    gps2m$angle.rad[i+1] <- atan2(b, a)
}
# convert to angle data
gps2m$angle.rad <- as.trig(gps2m$angle.rad)  # set it as angles
gps2m$bearing_gps <- as.bearing(gps2m$angle.rad)  # convert to 360° angles
gps2m <- rename(gps2m, lon_gps = lon, lat_gps = lat)


# Compute the mean bearing over 15s from the GPS GPGGA data
gps2m$group <- rep(1:round_any(nrow(gps2m)/15, 1), each = 15)[-nrow(gps2m)-1]




# Extract the same time period of the TS data
ts2m <- filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:09:50"))
ts2m$angle.rad_ts <- NA


# Compute bearing between each point and the next one (1s)
for (i in 1:(nrow(ts2m)-1)) {
    a <- geodDist(ts2m$lat[i+1], ts2m$lon[i+1], ts2m$lat[i], ts2m$lon[i])
    b <- geodDist(ts2m$lat[i], ts2m$lon[i], ts2m$lat[i+1], ts2m$lon[i])
    # check if positive values and correct if 
    if (ts2m$lon[i+1] < ts2m$lon[i]) {
    a <- -a
    }
    if (ts2m$lat[i+1] < ts2m$lat[i]) {
    b <- -b
    }
    ts2m$angle.rad_ts[i+1] <- atan2(b, a)
}
# convert to angle data
ts2m$angle.rad_ts <- as.trig(ts2m$angle.rad_ts)  # set it as angles
ts2m$bearing_ts <- as.bearing(ts2m$angle.rad_ts)  # convert to 360° angles




# Join the two traj with their bearing
ts_gps <- join(ts2m, gps2m)

select(ts_gps, bearing, bearing_gps)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TS RAW BEARINGS ARE UNUSABLE
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

select(ts_gps, bearing_ts, bearing_gps)
# But calculated bearing seems ok. 




# Join GPS with TS data

# Some means are biased by just one point (e.g. mean = 310°, and 1 is 40°)
# But it's very similar on average which is good. 
# We can add a filter that removes the data for which the diff in bearing is too high and due to a lack of accuracy of gps data


ts_gps <- cbind(select(ts2m, lon, lat, bearing_ts, dateTimeUTC), ddply(gps2m, ~group, function(x) { 
	# x <- gps2m[which(gps2m$group == 7), ]
	
	# FILTER BIASED CHANGES IN BEARING
	x <- x[which(diff(x$bearing_gps) < abs(30)), ]
	mean = mean.circular(as.bearing(x$bearing_gps))
	lon = max(x$lon_gps)
	lat = min(x$lat_gps)
	
	data.frame(mean_bearing_gps = mean, lon_gps = lon, lat_gps = lat)
}))


# Compute difference between mean bearing from gps and from ts
mean(na.omit(as.numeric(ts_gps$mean_bearing_gps) - as.numeric(ts_gps$bearing_ts)))
# MEAN OF -0.1147336° OF DIFFERENCE BETWEEN THE TOO. GOOD. 


# # Compare gps from TS and gps GPGGA visually
# ggplot() + geom_point(data = ts_gps_2m, aes(x = lon, y = lat, colour = bearing_ts)) + geom_point(data = ts_gps_2m, aes(x = lon_gps, y = lat_gps, colour = mean_bearing))

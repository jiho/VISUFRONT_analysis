#
#      Read and cut GPS data
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("oce")
library("dplyr")
library("ggplot2")
library("lubridate")
library("discuss")


# Concatenate all GPS files
gpsFiles <- list.files("/Users/faillettaz/Dropbox/visufront-data/gps/", pattern=glob2rx("gps*"), full=TRUE)
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
dateTimeUTC <- ymd_hms(dateTime)

# compensate for time difference with computer time
# dateTime <- dateTimeUTC + 2 * 3600 - 7


# Cleanup unused columns
gpsDf <- data.frame(lat, lon, dateTimeUTC)
head(gpsDf)
tail(gpsDf)

gpsDf <- filter(gpsDf, !is.na(lon))


# Cut into transects
transects <- read.csv("transects.csv", na.strings="", sep = ";", colClasses=c("character", "POSIXct", "POSIXct"))

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
filter(gpsLag, dateTimeUTC > as.POSIXct("2013-07-27 03:40:00") & dateTimeUTC < as.POSIXct("2013-07-27 09:10:00"))
# NO DATA BETWEEN 3:42 AM to 9:05 AM on the 27th
# 43.55036 7.586471 2013-07-27 03:42:55
# 43.60356 7.417247 2013-07-27 09:05:40


# Check very fine scale path 
ggplot(data = filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:00") & dateTimeUTC < ymd_hms("2013-07-26 22:01:00")), aes(x = lon, y = lat)) + geom_path() + geom_point()
# Very straight, this is good. 




gpsLag[which(gpsLag$dateTimeUTC > ymd_hms("2013-07-26 23:30:00") & gpsLag$dateTimeUTC < as.POSIXct("2013-07-26 23:31:00")), ]

filter(gpsLag, dateTimeUTC > as.POSIXct("2013-07-26 23:30:00"))


, dateTimeUTC < as.POSIXct("2013-07-26 23:31:00"))


filter(gpsLag, dateTimeUTC > as.POSIXct("2013-07-27 23:30:00"))[1:100,]


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
tsLag <- filter(ts, dateTimeUTC > as.POSIXct("2013-07-26 17:07:00") & dateTimeUTC < as.POSIXct("2013-07-28 02:10:00"))
# For some reasons, the time are wrong
ggplot() + geom_path(data= tsLag, aes(x = lon, y = lat))

# Check bearings distribution
ggplot() + geom_bar(aes(x = bearing), data = tsLag, binwidth = 5) + polar()


# Check very fine scale path 
ggplot(data= filter(tsLag, dateTimeUTC > as.POSIXct("2013-07-26 18:07:00") & dateTimeUTC < as.POSIXct("2013-07-26 18:09:00")), aes(x = lon, y = lat)) + geom_path() + geom_point()
# Very straight, this is good. 












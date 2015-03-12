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
library("gridExtra")
library("geosphere")


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


# # Cut into transects
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
#




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
tsFiles <- list.files(str_c(data_dir,"_raw_/TS"), pattern="*.tethys", full=TRUE)

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
tsLag <- filter(ts, dateTimeUTC > ymd_hms("2013-07-26 21:30:00") & dateTimeUTC < ymd_hms("2013-07-28 05:24:00"))
ggplot() + geom_path(data= tsLag, aes(x = lon, y = lat))

# Check bearings distribution
ggplot() + geom_bar(aes(x = bearing), data = tsLag, binwidth = 5) + polar()


# Check very fine scale path 
ggplot(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 21:33:00") & dateTimeUTC < ymd_hms("2013-07-26 21:35:00")), aes(x = lon, y = lat)) + geom_path() + geom_point()
# Very straight, this is good.


# Check path and bearing
ggplot(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 21:33:00") & dateTimeUTC < ymd_hms("2013-07-28 05:24:00")), aes(x = lon, y = lat)) + geom_path() + geom_point(aes(color = bearing)) + coord_quickmap()




# Compare gps from TS and gps GPGGA
ggplot() + 
# TS data
geom_path(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:03:50")), aes(x = lon, y = lat), colour = "red") + geom_point(data = filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:03:50")), aes(x = lon, y = lat), colour = "red") +
# GPS GPGGA data
geom_path(data = filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:03:50")), aes(x = lon, y = lat)) + geom_point(data = filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:03:50")), aes(x = lon, y = lat)) + coord_map() #+ coord_cartesian(xlim = c(7.56, 7.562), ylim = c(43.632, 43.634))
# ggsave(file = "plot/gps_TS_GPGAA.pdf")




# Compute bearing of from the GPS data for the 2 min period
# --------------------------------------------------------------------
gps2m <- filter(gpsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:03:50"))

# Compute bearing between each point and the next one (1s)
gps2m$bearing_gps <- NA
for (i in 1:(nrow(gps2m)-1)) {
	gps2m$bearing_gps[i] <- bearing(select(gps2m, lon, lat)[i,], select(gps2m, lon, lat)[i+1,])
}



# Extract the same time period of the TS data
ts2m <- filter(tsLag, dateTimeUTC > ymd_hms("2013-07-26 22:00:04") & dateTimeUTC < ymd_hms("2013-07-26 22:03:50"))

# Compute bearing between each point and the next one (1s)
ts2m$bearing_ts <- NA
for (i in 1:(nrow(ts2m)-1)) {
	ts2m$bearing_ts[i] <- bearing(select(ts2m, lon, lat)[i,], select(ts2m, lon, lat)[i+1,])
}



# # Compare visually
# ggplot() + geom_point(data = ts2m, aes(x = lon, y = lat, colour = bearing_ts)) + geom_point(data = gps2m, aes(x = lon, y = lat, colour = as.numeric(bearing_gps))) + coord_map()

# Compare means
# GPS FROM TS and RAW BEARING FROM TS
mean(na.omit(ts2m$bearing_ts)) - mean(na.omit(ts2m$bearing)) # ~3.7°
# GPS GPGGA and RAW BEARING FROM TS
mean(na.omit(gps2m$bearing_gps)) - mean(na.omit(ts2m$bearing)) # ~2.5°
# GPS GPGGA and GPS FROM TS
mean(na.omit(gps2m$bearing_gps)) - mean(na.omit(ts2m$bearing_ts)) # ~1.2°

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TS RAW BEARINGS ARE USABLE but less precise
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







# COMPUTE MEAN BEARING PER MIN AND CREATE TABLE FILE 
# --------------------------------------------------------------------

# Split by transect
transects <- read.csv("transects.csv", na.strings="", sep = ";", colClasses=c("character", "POSIXct", "POSIXct"))[15:28, ]


# Subset data to work on fewer first
gpsLag <- filter(gpsDf, dateTimeUTC > ymd_hms("2013-07-26 21:30:00") & dateTimeUTC < ymd_hms("2013-07-28 05:25:00"))
range(gpsLag$dateTimeUTC)



gps1mAll <- adply(transects$name, 1, function(x) {
	
	# x <- transects$name[1]
	transect <- filter(transects, name == x)
	
	gps <- filter(gpsLag, dateTimeUTC > (ymd_hms(transect$dateTimeStart) - 2 * 3600 - 60) & dateTimeUTC < (ymd_hms(transect$dateTimeEnd) - 2 * 3600 + 60))
	
	# if 1s gps data are not available, read gps from TS
	if (transect$name %in% c("lagrangian_3", "lagrangian_4", "lagrangian_5")) {
		gps <- filter(ts, dateTimeUTC > (ymd_hms(transect$dateTimeStart) - 2 * 3600 - 60) & dateTimeUTC < (ymd_hms(transect$dateTimeEnd) - 2 * 3600 + 60))
	}

	# Compute bearing between each point and the next one (1s)
	gps$bearing <- NA
	for (i in 1:(nrow(gps)-1)) {
		gps$bearing[i] <- bearing(select(gps, lon, lat)[i,], select(gps, lon, lat)[i+1,])
	}


	# Compute the mean bearing over 1m from the GPS GPGGA data
	gps$dateTimeUTC_min <- round_any(gps$dateTimeUTC, 60)
	
	# Remove the last line that will never have a bearing
	gps <- gps[-nrow(gps), ]
	
	gps_1m <- ddply(gps, ~dateTimeUTC_min, function(x) {
		
		# print(unique(x$dateTimeUTC_min))
		
		# x <- gps[which(gps$dateTimeUTC_min == ymd_hms("2013-07-26 21:33:00")), ]

		# filter large changes in bearing within minutes
		if (nrow(x) > 1 & any(diff(x$bearing_gps) < abs(30))) {
	  	x <- x[which(diff(x$bearing_gps) < abs(30)), ]
		}
		
		# Compute circular mean and sd
	  	mean <- mean.circular(as.bearing(x$bearing))
		sd <- sd.circular(as.bearing(x$bearing))
		
		
		# Get the mean lon and lat instead of the min of max no to loose the beginning and end
		lon <- mean(x$lon)
		lat <- mean(x$lat)
	
	  	data.frame(lon = lon, lat = lat, mean_bearing_gps = mean, sd = sd)
		
	  })

}, .progress = "text")
gps1mAll <- gps1mAll[, -1]
head(gps1mAll)
# tail(gps1mAll)
# gps1mAll[1:100,]




# Interpolate bearing and lon lat to get 1s time steps
# --------------------------------------------------------------------

# Interpolate a slice of data for which the x-axis is time
interp.time <- function(x, y, y.name, x.step=1,  ...) {

  # interpolate as numbers
  time <- as.numeric(x)
  i <- approx(x = time, y = y, xout = seq(min(time), max(time), by=x.step))

  # extract a data.frame
  out <- data.frame(dateTimeSec = i$x, var = i$y)
  colnames(out) <- c("dateTimeSec", y.name)

  # reconvert to time
  out$dateTimeSec <- as.POSIXct(out$dateTimeSec, origin=as.POSIXct("1970-01-01 01:00:00.00"))
  # range(out$dateTimeSec)
  return(out)
}



# Run interpolation of lon, lat, bearing and sd per second
gpsAll <- adply(transects$name[c(3,4,5)], 1, function(x) {
	
	# transect <- transects[2,]
	transect <- filter(transects, name == x)
	# print(transect)
	
	# Select data (-+60 s to start having good data at the begining of the transects)
	gps <- gps1mAll[which(gps1mAll$dateTimeUTC_min > (ymd_hms(transect$dateTimeStart) - 2 * 3600 - 60) & gps1mAll$dateTimeUTC_min < (ymd_hms(transect$dateTimeEnd) - 2 * 3600 + 30)), ]
	
	# Interpolation
	lontmp <- interp.time(x = gps$dateTimeUTC_min, y = gps$lon, x.step = 1, y.name = "lon")
	lattmp <- interp.time(x = gps$dateTimeUTC_min, y = gps$lat, x.step = 1, y.name = "lat")
	bearingtmp <- interp.time(x = gps$dateTimeUTC_min, y = gps$mean_bearing_gps, x.step = 1, y.name = "bearing")
	bearingsdtmp <- interp.time(x = gps$dateTimeUTC_min, y = gps$sd, x.step = 1, y.name = "sd")

	# format data
	data.frame(transect = rep(transect$name, nrow(lontmp)), dateTimeUTC = lontmp$dateTimeSec, lon = round_any(lontmp$lon, 0.00001), lat = round_any(lattmp$lat, 0.00001), bearing = round_any(bearingtmp$bearing, 0.001), sd = round_any(bearingsdtmp$sd, 0.005))
})
gpsAll <- gpsAll[, -1]
head(gpsAll)
# tail(gpsAll)



# inspect visually
ggplot() + geom_point(data= gpsAll, aes(x = lon, y = lat, colour = bearing)) + geom_path(data = tsLag, aes(x = lon, y = lat)) + coord_map()
# ggsave("traj_bearing.pdf")



ddply(gpsAll, ~transect, summarize, mean = mean(bearing), sd = mean(sd))
#         transect      mean         sd
# 1   lagrangian_1 331.88962 0.08937468
# 2   lagrangian_2 153.10690 0.07772209
# 3   lagrangian_3 333.00891 0.01130073
# 4   lagrangian_4 171.87109 0.01670370
# 5   lagrangian_5 342.12104 0.01067944
# 6   lagrangian_6 162.39459 0.07202128
# 7   lagrangian_7 338.51732 0.09358875
# 8   lagrangian_8 161.56685 0.11271465
# 9   lagrangian_9 331.38366 0.09541660
# 10 lagrangian_11 180.59167 0.07172279
# 11 lagrangian_12  11.91548 0.07700975
# 12 lagrangian_13 202.45207 0.08436959
# 13 lagrangian_14  20.74815 0.07068693
# 14 lagrangian_15 205.71080 0.06690580






# Write TXT file
# --------------------------------------------------------------------

# bearing per seconde
write.table(select(gpsAll, -transect), file = "bearing-visufront-from-GPS-1SEC.txt", fileEncoding = "ASCII", append = F, row.names = F, sep = ";")




# bearing per minute
gps1mAll$lon <- round_any(gps1mAll$lon, 0.00001)
gps1mAll$lat <- round_any(gps1mAll$lat, 0.00001)
gps1mAll$bearing <- as.numeric(round_any(gps1mAll$mean_bearing_gps, 0.001))
gps1mAll$sd <- round_any(gps1mAll$sd, 0.001)
gps1mAll <- select(rename(select(gps1mAll, -mean_bearing_gps), dateTimeUTC = dateTimeUTC_min), dateTimeUTC, lon, lat, bearing, sd)

write.table(gps1mAll, file = "bearing-visufront-from-GPS-1MIN.txt", fileEncoding = "ASCII", append = F, row.names = F, sep = ";")





# Extract raw bearing from the TS GPS (every 15s)
tsLag15s <- adply(transects$name, 1, function(x) {
	
	# transect <- transects[2,]
	transect <- filter(transects, name == x)
	# print(transect)
	
	# Select data (-+60 s to start having good data at the begining of the transects)
	gps <- tsLag[which(tsLag$dateTimeUTC > (ymd_hms(transect$dateTimeStart) - 2 * 3600 - 60) & tsLag$dateTimeUTC < (ymd_hms(transect$dateTimeEnd) - 2 * 3600 + 30)), ]
	
	# Compute bearing between each point and the next one (1s)
	gps$bearing <- NA
	for (i in 1:(nrow(gps)-1)) {
		gps$bearing[i] <- round_any(bearing(select(gps, lon, lat)[i,], select(gps, lon, lat)[i+1,]), 0.001)
	}
	select(gps, dateTimeUTC, lon, lat, bearing)
})
head(tsLag15s)
tail(tsLag15s)

# Format data
tsLag15s <- tsLag15s[-nrow(tsLag15s), -1]
tsLag15s <- cbind(select(tsLag15s, dateTimeUTC, lon, lat, bearing), sd = rep(0, nrow(tsLag15s)))



write.table(tsLag15s, file = "bearing-visufront-from-TS-15SEC.txt", fileEncoding = "ASCII", append = F, row.names = F, sep = ";")






# Extract raw position data and interpolate per 1s FOR TS
# --------------------------------------------------------------------
Ts1s <- adply(transects$name[c(3,4,5)], 1, function(x) {
	
	# transect <- transects[3,]
	transect <- filter(transects, name == x)
	# print(transect)
	
	# Select data (-+60 s to start having good data at the begining of the transects)
	gps <- tsLag[which(tsLag$dateTimeUTC > (ymd_hms(transect$dateTimeStart) - 2 * 3600 - 60) & tsLag$dateTimeUTC < (ymd_hms(transect$dateTimeEnd) - 2 * 3600 + 30)), ]
	
	# Compute bearing between each point and the next one (1s)
	gps$bearing_gps <- NA
	for (i in 1:(nrow(gps)-1)) {
		gps$bearing_gps[i] <- bearing(select(gps, lon, lat)[i,], select(gps, lon, lat)[i+1,])
	}
	
	# add true bearing if available
	if (!is.na(gps$bearingTrue[nrow(gps)]) & is.na(gps$bearing[nrow(gps)])) {
	gps$bearing_gps[nrow(gps)] <- gps$bearingTrue[nrow(gps)]
	}
	
	if (is.na(gps$bearingTrue[nrow(gps)]) & !is.na(gps$bearing[nrow(gps)])) {
	gps$bearing_gps[nrow(gps)] <- gps$bearing[nrow(gps)]
	}
	
	
	# Interpolation
	lontmp <- interp.time(x = gps$dateTimeUTC, y = gps$lon, x.step = 1, y.name = "lon")
	lattmp <- interp.time(x = gps$dateTimeUTC, y = gps$lat, x.step = 1, y.name = "lat")
	bearingtmp <- interp.time(x = gps$dateTimeUTC, y = gps$bearing_gps, x.step = 1, y.name = "bearing")

	# format data
	data.frame(transect = rep(transect$name, nrow(lontmp)), dateTimeUTC = lontmp$dateTimeSec, lon = round_any(lontmp$lon, 0.00001), lat = round_any(lattmp$lat, 0.00001), bearing = round_any(bearingtmp$bearing, 0.001))
})
Ts1s <- Ts1s[, -1]
head(Ts1s)
tail(Ts1s)

Ts1s <- filter(Ts1s, !is.na(bearing))



# Extract raw position data per 1s WITH GPS
# --------------------------------------------------------------------
gps1s <- adply(transects$name[-c(3,4,5)], 1, function(x) {
	
	transect <- transects[2,]
	transect <- filter(transects, name == x)
	# print(transect)
	
	# Select data (-+60 s to start having good data at the begining of the transects)
	gps <- gpsLag[which(gpsLag$dateTimeUTC > (ymd_hms(transect$dateTimeStart) - 2 * 3600 - 60) & gpsLag$dateTimeUTC < (ymd_hms(transect$dateTimeEnd) - 2 * 3600 + 30)), ]
	
	# Compute bearing between each point and the next one (1s)
	gps$bearing_gps <- NA
	for (i in 1:(nrow(gps)-1)) {
		gps$bearing_gps[i] <- bearing(select(gps, lon, lat)[i,], select(gps, lon, lat)[i+1,])
	}
		

	# format data
	data.frame(transect = rep(transect$name, nrow(gps)), dateTimeUTC = gps$dateTimeUTC, lon = round_any(gps$lon, 0.00001), lat = round_any(gps$lat, 0.00001), bearing = round_any(gps$bearing_gps, 0.001))
}, .progress = "text")
gps1s <- gps1s[, -1]
head(gps1s)
tail(gps1s)
gps1s <- filter(gps1s, !is.na(bearing))




# Combine raw GPS 1s with interp TS 1s
# --------------------------------------------------------------------
bearings <- filter(arrange(rbind(gps1s, Ts1s), dateTimeUTC), !is.na(bearing))

# bearing per seconde
write.table(select(bearings, -transect), file = "bearing-visufront-GPS_RAW_1SEC-TS_INTERP_1SEC.txt", fileEncoding = "ASCII", append = F, row.names = F, sep = ";")



# Combine raw GPS 1s and raw TS 15s
# --------------------------------------------------------------------
bearingsRaw <- filter(arrange(rbind(select(gps1s, -transect), tsLag15s), dateTimeUTC), !is.na(bearing))
head(bearingsRaw)
write.table(bearingsRaw, file = "bearing-visufront-GPS_RAW_1SEC-TS_RAW_15SEC.txt", fileEncoding = "ASCII", append = F, row.names = F, sep = ";")




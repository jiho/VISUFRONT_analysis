#
#      Functions to process raw data
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

# Read thermosalinometer from Tethys
read.ts <- function(file) {
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
  midnight <- as.POSIXct(str_c(year, "-", month, "-", day, " 00:00:01"))
  d$dateTime <- midnight + d$timeUTC
  # convert into local time
  d$dateTime <- d$dateTime + 2*3600  
  
  # fill the gaps in GPS coordinates when necessary
  if (any(is.na(d$lat))) {
    d$lat <- approx(x=d$dateTime, y=d$lat, xo=d$dateTime)$y
    d$lon <- approx(x=d$dateTime, y=d$lon, xo=d$dateTime)$y
  }
  
  # select interesting data columns
  d <- d[,c("dateTime", "lat", "lon", "atmPressure", "tempAir", "humidity", "windDirection", "windSpeed", "depth", "temperature", "salinity", "fluorometry")]
	
	return(d)
}

# Read hydrological data from ISIIS
read.isiis <- function(file) {
	options(digits.secs=2)  # allow split seconds
	
	# read the data
  d <- read.delim(file, skip=10, fileEncoding="ISO-8859-1", encoding="UTF-8", stringsAsFactors=FALSE)

  # clean names
	# remove double dots
  names(d) <- str_replace_all(names(d), fixed(".."), ".")
  names(d) <- str_replace_all(names(d), fixed(".."), ".")
	# remove dots at end of names
  names(d) <- str_replace(names(d), "\\.$", "")

  # extract date from file name
  year <- str_sub(file, -16, -13)
  month <- str_sub(file, -12, -11)
  day <- str_sub(file, -10, -9)

  # compute date and time
  d$dateTimeMsec <- as.POSIXct(str_c(year, "-", month, "-", day, " ", d$Time))
  # detect midnight shifts
  midnightShift <- which(diff(d$dateTimeMsec) < 0)
  if (length(midnightShift) > 0) {
    d$dateTimeMsec[midnightShift:nrow(d)] <- d$dateTimeMsec[midnightShift:nrow(d)] + 24 * 3600
  }

  # keep important columns
  d <- d[,c("dateTimeMsec", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm")]
	
	return(d)
}
dist.from.shore <- function(lat, lon) {
	# find the most northern-wastern point
	pointLat <- max(lat, na.rm=TRUE)
	pointLon <- min(lon, na.rm=TRUE)
	
	# compute distance from shore for all points
	coast <- read.csv("map/gshhg_coteazur_i.csv")
	
	# find the closest point on the coast
	library("oce")
	dists <- geodDist(lat1=coast$lat, lon1=coast$lon, lat2=pointLat, lon2=pointLon)
	minIdx <- which.min(dists)
	refLat <- coast$lat[minIdx]
	refLon <- coast$lon[minIdx]
	
	# now compute distance
	dists <- geodDist(lat1=lat, lon1=lon, lat2=refLat, lon2=refLon)
	
	return(dists)
}







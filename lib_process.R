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
  library("stringr")
  library("lubridate")

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
  d$dateTimeMsec <- ymd_hms(str_c(year, "-", month, "-", day, " ", d$Time))
  # detect midnight shifts
  midnightShift <- which(diff(d$dateTimeMsec) < 0)
  if (length(midnightShift) > 0) {
    d$dateTimeMsec[midnightShift:nrow(d)] <- d$dateTimeMsec[midnightShift:nrow(d)] + 24 * 3600
  }

  # keep important columns
  d <- d[,c("dateTimeMsec", "Pressure.dbar", "Depth.m", "Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm")]

	return(d)
}

# Detect up and down casts in a depth yo
detect.casts <- function(depth, order=200) {
  # smoothing the depth profile using a moving average and find the turning points

	# smooth depths
  library("pastecs")
  depth_avg <- decaverage(-depth, times=3, weights=c(seq(1, order), order+1, seq(order, 1, -1)))
  # plot(depth_avg)
  depth_avg <- as.numeric(pastecs::extract(depth_avg, component="filtered"))

  # detect turning points
  TP <- suppressWarnings(turnpoints(depth_avg))

  # set cast numbers (different for up and down casts)
  cast <- cumsum(TP$peaks | TP$pits) + 1

  # detect which are up and which are down casts:
  # if the first turning point is a peak, then the first cast (and all odd casts) are upcasts
  if ( TP$firstispeak ) {
    # these are the types for
    #              even  & odd   cast numbers
    castTypes <- c("down", "up")
  } else {
    castTypes <- c("up", "down")
  }
  down.up <- castTypes[cast %% 2 + 1]

	return(data.frame(cast, down.up))
}

# Compute the straight line distance from the starting point of a lat,lon trajectory
dist.from.start <- function(lat, lon) {
	library("oce")
	geodDist(lon1=lon, lat1=lat, lon2=na.omit(lon)[1], lat2=na.omit(lat)[1])
}

# Compute the distance from Villefranche
dist.from.villefranche <- function(lat, lon) {
	# TODO should compensate for the length of cable put out and the angle of the cable (i.e. depth of ISIIS)
	library("oce")
	geodDist(lon1=lon, lat1=lat, lon2=7.3118057, lat2=43.70528)
}

# Compute the distance from shore, for the VISUFRONT cruise only
dist.from.shore <- function(lat, lon) {
	# find the most north-western point of the transect (the one closest to shore already)
	base_point_lat <- max(lat, na.rm=TRUE)
	base_point_lon <- min(lon, na.rm=TRUE)

	# find the point on the coast closest to this point
	# TODO should find the point of intersection between the coast and the axis of the transect
	coast <- read.csv("map/gshhg_coteazur_i.csv")
	library("oce")
	dists <- geodDist(lon1=coast$lon, lat1=coast$lat, lon2=base_point_lon, lat2=base_point_lat)
	minIdx <- which.min(dists)
	ref_lat <- coast$lat[minIdx]
	ref_lon <- coast$lon[minIdx]

	# compute distances from this point to every point in the track
	# TODO should compensate for the length of cable put out and the angle of the cable (i.e. depth of ISIIS)
	dists <- geodDist(lon1=lon, lat1=lat, lon2=ref_lon, lat2=ref_lat)

	return(dists)
}

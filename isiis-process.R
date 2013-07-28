#
#      Read and process ISIIS physical data
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
library("oce")

source("lib_process.R")	# dist.from.shore

# get data
hydroFiles <- list.files("ISIIShydro", pattern=glob2rx("ISIIS*.txt"), full=T)

# read data from all transects
options(digits.secs=2)  # allow split seconds
d <- adply(hydroFiles, 1, function(file) {
  d <- read.delim(file, skip=10, fileEncoding="ISO-8859-1", encoding="UTF-8", stringsAsFactors=FALSE)

  # clean names
  names(d) <- str_replace_all(names(d), fixed(".."), ".")
  names(d) <- str_replace_all(names(d), fixed(".."), ".")
  # names(d) <- str_replace(names(d), perl("\\\\.$"), ".") # does not work

  # extract date from file name
  date <- str_replace(file, "ISIIShydro/ISIIS", "")
  date <- str_replace(date, ".txt", "")
  year <- str_sub(date, 1, 4)
  month <- str_sub(date, 5, 6)
  day <- str_sub(date, 7, 8)

  # compute date and time
  d$dateTimeMsec <- as.POSIXct(str_c(year, "-", month, "-", day, " ", d$Time))
  # detect midnight shifts
  midnightShift <- which(diff(d$dateTimeMsec) < 0)
  if (length(midnightShift) > 0) {
    d$dateTimeMsec[midnightShift:nrow(d)] <- d$dateTimeMsec[midnightShift:nrow(d)] + 24 * 3600
  }

  # keep important columns
  d <- d[,c("dateTimeMsec", "Pressure.dbar.", "Depth.m.", "Temp.C.", "Salinity.PPT.", "Fluoro.volts.", "Oxygen.ml.l.", "Irrandiance.UE.cm.")]

  return(d)
}, .progress="text")
d <- d[,-1]


##{ Cleanup data ----------------------------------------------------------

# remove impossible values
d$Temp.C.[d$Temp.C. <= 0] <- NA
d$Salinity.PPT.[d$Salinity.PPT. <= 10] <- NA
d$Fluoro.volts.[d$Fluoro.volts. <= 0.01] <- NA
d$Oxygen.ml.l.[d$Oxygen.ml.l. <= 0] <- NA

# check TS diagram
ggplot(d) + geom_point(aes(x=Temp.C., y=Salinity.PPT.), size=1, alpha=0.1, na.rm=T)

# check all variables
dm <- melt(d, id.vars=c("dateTimeMsec", "Pressure.dbar.", "Depth.m."))
ggplot(dm) + geom_point(aes(x=dateTimeMsec, y=value), size=1, alpha=0.1, na.rm=T) + facet_wrap(~variable, scale="free_y")

# add GPS data
ts <- read.csv("ts.csv", stringsAsFactors=FALSE)
ts$dateTime <- as.POSIXct(ts$dateTime)
# base ourselves on the second, it is enough already
d$dateTime <- round(d$dateTimeMsec)
d <- join(d, ts[,c("dateTime", "lat", "lon")], by="dateTime")

# interpolate GPS data
sum(is.na(d$lat))
sum(is.na(d$lon))
d$lat <- approx(x=as.numeric(d$dateTime), y=d$lat, xo=as.numeric(d$dateTime))$y
d$lon <- approx(x=as.numeric(d$dateTime), y=d$lon, xo=as.numeric(d$dateTime))$y

# using salinity, temperature and pressure for calculation of seawater density
# using UNESCO formulation
d$Density <- swRho(d$Salinity.PPT., d$Temp.C., d$Pressure.dbar., eos="unesco")

# }


##{ Cut by transect -------------------------------------------------------

isiis <- d
write.csv(isiis, file="isiis.csv", row.names=FALSE)

transects <- read.csv("transects.csv", na.strings=c("", "NA"), colClasses=c("character", "POSIXct", "POSIXct"))

pdf("isiis-transects.pdf")
d_ply(transects, ~name, function(x, data) {
	message(x$name)

  # extract the portion of the data corresponding to this transect
  cData <- data[which(data$dateTime > x$dateTimeStart-5 & data$dateTime < x$dateTimeEnd+5),]

  if (nrow(cData) >= 1) {
    
    # compute distance from first point and from a reference point
    cData$distanceFromStart <- geodDist(lat1=cData$lat, lon1=cData$lon, lat2=na.omit(cData$lat)[1], lon2=na.omit(cData$lon)[1]) / 1.852
    cData$distanceFromVlfr <- geodDist(cData$lat, cData$lon, lat2=43.70528, lon2=7.3118057) / 1.852
    cData$distanceFromShore <- dist.from.shore(cData$lat, cData$lon) / 1.852
    
    # detect up and down casts by smoothing the depth profile and finding the turning points
    # smooth depths using a moving average
    library("pastecs")
    order <- 200
    depth_avg <- decaverage(-cData$Depth.m., times=3, weights=c(seq(1, order), order+1, seq(order, 1, -1)))
    # plot(depth_avg)
    depth_avg <- as.numeric(pastecs::extract(depth_avg, component="filtered"))

    # detect turning points
    TP <- suppressWarnings(turnpoints(depth_avg))
    
    # give cast numbers (different for up and down casts)
    cData$cast <- cumsum(TP$peaks | TP$pits) + 1
    
    # detect which are up and which are down casts
    # if the first turning point is a peak, then the first cast (and all odd casts) are upcasts
    if ( TP$firstispeak ) {
      # these are the types for
      #              even  & odd   cast numbers
      castTypes <- c("down", "up")
    } else {
      castTypes <- c("up", "down")
    }
    cData$down.up <- castTypes[cData$cast %% 2 + 1]
    
    # plot to check
    print(ggplot(cData) + geom_point(aes(x=distanceFromShore, y=-Depth.m., colour=down.up), na.rm=T) + ggtitle(x$name))
  
    # store data file
    dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
    dataName <- deparse(substitute(data))
    write.csv(cData, file=str_c("transects/", x$name, "/", dataName, ".csv"), row.names=FALSE)    
  }
  
  return(NA)

}, data=isiis)
dev.off()

# }

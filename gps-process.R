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


# Concatenate all GPS files
gpsFiles <- list.files("ISIIShydro", pattern=glob2rx("gps*"), full=TRUE)
gps <- alply(gpsFiles, 1, function(file) {
  # read the content of the file
  content <- scan(file, what="character", sep="\n", fileEncoding="ISO-8859-1", encoding="UTF-8")
  return(content)
})
gps <- unlist(gps)

# only keep GPZDA (date) and GPGLL (lat, lon) lines
gpsD <- gps[str_detect(gps, fixed("$GPZDA"))]
gpsL <- gps[str_detect(gps, fixed("$GPGLL"))]

# # fix problems
# gpsD <- str_replace(gpsD, fixed("07<2013"), "07,2013")
# gpsD <- str_replace(gpsD, fixed("+0\u0082àR¢j"), "+00,00")
# gpsD <- str_replace(gpsD, fixed("$GPGGA,072104.00,4325.99442,N,00747.22920,E,1,09,0.9,17.39,M"), "")


# Read concatenated data
read.gps <- function(x) {
  temp <- tempfile()
  cat(x, file=temp, sep="\n")
  x <- read.table(temp, sep=",", header=FALSE, stringsAsFactors=FALSE, colClasses="character")
  file.remove(temp)
  return(x)
}

# gpsD <- read.gps(gpsD)
gpsL <- read.gps(gpsL)


# Extract coordinates
# the current format is DDMM.MMM
deglat <- as.numeric(str_sub(gpsL$V2, 1, 2))
minlat <- as.numeric(str_sub(gpsL$V2, 3))
deglon <- as.numeric(str_sub(gpsL$V4, 1, 3))
minlon <- as.numeric(str_sub(gpsL$V4, 4))
lat <- deglat + minlat/60
lon <- deglon + minlon/60

# Extract time
day <- str_sub(gpsL$V1, 10, 11)
month <- "07"
year <- "2013"

hour <- str_sub(gpsL$V6, 1, 2)
min <- str_sub(gpsL$V6, 3, 4)
sec <- str_sub(gpsL$V6, 5, 6)

dateTime <- str_c(str_c(year, month, day, sep="-"), " ", str_c(hour, min, sec, sep=":"))
dateTime <- as.POSIXct(dateTime)
# compensate for time difference with computer time
dateTime <- dateTime + 2 * 3600 - 7


# Cleanup unused columns
gps <- data.frame(lat, lon, dateTime)

# Cut into transects
transects <- read.csv("transects.csv", na.strings="", colClasses=c("character", "POSIXct", "POSIXct"))

d_ply(transects, ~name, function(x, gps) {
  # extract the portion of the GPS track
  cGps <- gps[gps$dateTime > x$dateTimeStart-5 & gps$dateTime < x$dateTimeEnd+5,]

  if (nrow(cGps) >= 1) {
    # compute distance from first point and from a reference point
    cGps$distanceFromStart <- geodDist(lat1=cGps$lat, lon1=cGps$lon, lat2=cGps$lat[1], lon2=cGps$lon[1]) * 1.852
    cGps$distanceFromVlfr <- geodDist(cGps$lat, cGps$lon, lat2=43.70528, lon2=7.3118057) * 1.852
  
    # store data file
    dir.create(str_c("transects/", x$name), showWarnings=FALSE, recursive=TRUE)
    write.csv(cGps, file=str_c("transects/", x$name, "/gps.csv"), row.names=FALSE)    
  }
  
}, gps=gps)

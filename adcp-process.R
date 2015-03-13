#
#     Fix bearings in ADCP record for VISUFRONT 
#
# (c) Copyright 2014 Jean-Olivier Irisson, GNU General Public License v3
#
#--------------------------------------------------------------------------

source("lib_process.R")
data_dir <- data_dir_path()

library("stringr")
library("plyr")
library("dplyr")


##{ Read all N1R files -----------------------------------------------------

# list N1R files
n1r_files <- list.files(str_c(data_dir, "/_raw_/ADCP172_000000.LTA/"), pattern=glob2rx("*.N1R"), full=TRUE)

# read all files
n1r <- ldply(n1r_files, function(f) {
  x <- scan(f, what="character", sep="\n", quiet=TRUE)
  data.frame(text=x, file=basename(f))
}, .progress="text")

# }


##{ Compute positions ------------------------------------------------------

# find lines with usable (= uninterrupted by ADCP) lat lon
xy_lno <- which(str_detect(n1r$text, "^\\$GPGLL") & !str_detect(n1r$text, "\\$PADCP"))

# parse all lat/lon
xy_text <- n1r$text[xy_lno]
xy <- read.table(text=paste(xy_text, collapse="\n"), sep=",")
# convert lat/lon from the GPGLL record into decimal lat/lon
decimal_lat_lon <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  deg <- round(x/100)
  min <- x - deg*100
  return(deg + min/60)
}
xy$lat <- decimal_lat_lon(xy[,2])
xy$lon <- decimal_lat_lon(xy[,4])
xy$lno <- xy_lno

# inspect positions
summary(xy)
# -> a few NAs

# detect how many lines are between each lat/lon record
table(diff(xy$lno))
# -> OK, not too far from each other

# simple plot (of a subsample of data)
library("ggplot2")
ggplot(xy[seq(1, nrow(xy), length=10000),]) + geom_point(aes(x=lon, y=lat), size=0.5)

# remove path outside the mission
xy_ok <- filter(xy, lon > 7)

summary(xy_ok)
# -> no more NAs! All were outside our region of interest

# }


##{ Compute missing bearings from positions --------------------------------

# find lines with missing bearings (either nothing or error code)
na_bearing_lno <- which(
  str_detect(n1r$text, "^\\$GPHDT,(999|666)\\.000") |
  str_detect(n1r$text, "^\\$GPHDT,-999.990") |
  str_detect(n1r$text, "^\\$GPHDT,,")
)

# for each missing bearing, find the previous line with a valid, non-NA lat-lon
index_of_previous <- function(x, y, ...) {
  # for each element of y, find the closest element of x which is before it and return its index
  floor(approx(x=x, y=1:length(x), xo=y, ...)$y)
}
# index_of_previous(x=c(10,20,30), y=c(22,34,5))
# index_of_previous(x=c(10,20,30), y=c(22,34,5), rule=2)
# index_of_previous(x=c(10,20,30), y=c(22,34,5), rule=c(1,2))

xy_i <- index_of_previous(x=xy$lno, y=na_bearing_lno)
any(is.na(xy_i))
# -> no NA = OK

# get the lines of the GPS locations before and after each missing bearing
xy_lno_before <- xy$lno[xy_i]
xy_lno_after  <- xy$lno[xy_i + 1]
# check that those lines are not too far from each other (to avoid computing bearings over large distances)
table(xy_lno_after - xy_lno_before)
# -> OK, most are not spread by very much

# get coordinates of points before and after
before_points <- xy[xy_i, c("lat", "lon")]
after_points <- xy[xy_i + 1, c("lat", "lon")]

# detect NA positions (the bearings function does not like them)
points <- cbind(before_points, after_points)
na_points <- is.na(rowSums(points))
points_no_na <- points[!na_points,]

# compute bearings between successive positions
library("geosphere")
bearings <- bearing(points_no_na[,1:2], points_no_na[,3:4])

# associate those bearings to the original line number
# NB: this recreates the NAs which is what we want
b <- data.frame(lno=na_bearing_lno)
b$bearing[!na_points] <- bearings

# at this point, the column "bearing" has some NAs:
# - some are created by NA longitudes
# - others are created by the fact that the before and after points are exactly at the same place

# replace missing value code by computed values in original text
# if computed bearing is NA, replace by a stardard missing value code
b$bearing[is.na(b$bearing)] <- 999

# replace missing value codes by computed bearings
text <- n1r$text[b$lno]
# str_replace(text, "\\$GPHDT,(999\\.000)*(?=(,|$)))", sprintf("$GPHDT,%07.3f,", b$bearing))
# replace completely unusual missing values with the usual missing value code, to ensure everything has the same width
text <- str_replace(text, fixed("$GPHDT,,"), "$GPHDT,999.000,")
text <- str_replace(text, fixed("$GPHDT,-999.990,"), "$GPHDT,999.000,")
# text <- str_replace(text, fixed("$GPHDT,666.000,"), "$GPHDT,999.000,") # NB: useless since 666.000 and 999.000 have the same width
# replace the missing value codes by the computed bearings
str_sub(text, 8, 14) <- sprintf("%07.3f", b$bearing)

n1r_ok <- n1r
text -> n1r_ok$text[b$lno]

# }


##{ Extract series of bearings ---------------------------------------------

# find lines with usable bearings
bearing_lno <- which(str_detect(n1r_ok$text, "^\\$GPHDT") & !str_detect(n1r_ok$text, "\\$PADCP"))

# parse all lat/lon
bearing_text <- n1r_ok$text[bearing_lno]
b <- read.table(text=paste(bearing_text, collapse="\n"), sep=",", na.strings=c("999.000"))
b <- b$V2

summary(b)
range(b, na.rm=T)

plot(bearings, type="l")
hist(bearings)

# looks OK
# TODO check for abrupt variation by a moving circular range window

# }


##{ Write the modified files -----------------------------------------------
  
d_ply(n1r_ok, ~file, function(x) {
  filename <- basename(x$file[1])
  cat(x$text, sep="\n", file=filename)
}, .progress="text")  

# }

#
#      Check transects delimitation
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#------------------------------------------------------------

f <- list.files("transects", recursive=T, pattern="isiis", full=T)

d <- alply(f, 1, read.csv, stringsAsFactors=FALSE)
transectName <- str_replace(f, "transects/", "")
transectName <- str_replace(transectName, "/isiis.csv", "")
attributes(d) <- NULL
names(d) <- transectName

ldply(d, function(x) {range(x$lat, na.rm=T)})
ldply(d, function(x) {range(x$lon, na.rm=T)})
ldply(d, function(x) {range(x$dateTime, na.rm=T)})

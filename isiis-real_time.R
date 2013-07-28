#
#      Read data from the ISIIS computer and Tethys on the fly and plot it
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")

source("lib_process.R")
source("lib_plot.R")

real.time.plot <- function(vars=c("Temp.C", "Salinity.PPT", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm")) {
	# get data
	hydroFiles <- list.files("/Volumes/ISIIShydro", pattern=glob2rx("ISIIS20130728*.txt"), full=T)
	d <- read.isiis(hydroFiles[length(hydroFiles)])

	tsFile <- "/Volumes/donnees/20130728.tethys"
	t <- read.ts(tsFile)

	# get position data
	d$dateTime <- round(d$dateTimeMsec)
	d <- join(d, t[,c("dateTime", "lat", "lon")], by="dateTime")
	# interpolate data points
	d$lat <- approx(x=as.numeric(d$dateTime), y=d$lat, xo=as.numeric(d$dateTime))$y
	d$lon <- approx(x=as.numeric(d$dateTime), y=d$lon, xo=as.numeric(d$dateTime))$y

	# compute distance from a reference point
	d$distanceFromVlfr <- dist.from.villefranche(d$lat, d$lon)

	# detect yos
	casts <- detect.casts(d$Depth.m)
	d <- cbind(d, casts)
	# check
	# ggplot(d) + geom_point(aes(x=dateTime, y=-Depth.m., colour=down.up))

	# interpolate all variables
	dm <- melt(d, id.vars=c("Depth.m", "distanceFromVlfr", "down.up"), measure.vars=vars)

	di <- ddply(dm, ~variable, function(x) {
	  x <- na.omit(x[which(x$down.up=="up"),])
	  xi <- interp.dist(x=x$distanceFromVlfr, y=x$Depth.m, z=x$value, duplicate="mean", x.step=300, y.step=1, anisotropy=1300)
	})
	di <- rename(di, c("x"="distance", "y"="Depth.m"))

	plots <- dlply(di, ~variable, function(x) {
	  ggplot(x, aes(x=distance, y=-Depth.m)) +
	    # geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
	    geom_tile(aes(fill=value), na.rm=T) +
	    stat_contour(aes(z=value), colour="white", alpha=0.7, bins=5, size=0.2, na.rm=TRUE) +
	    scale_fill_gradientn(colours=spectral(), guide="none", na.value=NA) +
	    scale_x_continuous(expand=c(0,0)) +
	    scale_y_continuous(expand=c(0,0))
	})

	do.call(grid.arrange, c(plots,list(ncol=1)))
}

while ( 1 == 1 ) {
	real.time.plot(c("Salinity.PPT", "Fluoro.volts"))
	Sys.sleep(60)
}

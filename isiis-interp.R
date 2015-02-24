#
#      Interpolate environmental data
#
#  (c) Copyright 2014 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

# TODO find good anisotropy ratio
# TODO read data from all source (ISIIS, CTD, gliders)
# TODO join them in a db by transect
# TODO interpolate all variables on all transects
# TODO - with a coarse grid for joining with other data
# TODO - with a fine grid (and possibly extrapolation) for plots

library("reshape2")
library("ggplot2")
library("plyr")
library("grid")
library("dplyr")
library("stringr")
library("gridExtra")
library("scales")

source("lib_plot.R")

# get all isiis data
isiisFiles <- list.files("transects", pattern="isiis.csv", full=TRUE, recursive=TRUE)
# isiisFiles <- "transects/cross_current_4/isiis.csv"

# interpolate all variables for each file
l_ply(isiisFiles, function(file) {

  # read data
  dir <- dirname(file)
  message(dir)
  e <- read.csv(file)
  
  # decide on which distance measure to use
  if ( str_detect(dir, "along") |  str_detect(dir, "lagrangian") ) {
    distance <- "dist_from_start"
  } else {
    distance <- "dist_from_shore"
  }

  # melt data to interpolate each variable sequentially
  em <- melt(e, id.vars=c(distance, "Depth.m", "cast", "down.up"), measure.vars=c("Salinity.PPT", "Temp.C", "Fluoro.volts", "Oxygen.ml.l", "Density"))
  # interp does not like NAs
  em <- na.omit(em)
  # homogenise distance name
  names(em)[1] <- "Distance.km"
  
  # bin over depth to avoid small scale variations which would throw off the interpolation
  em$Depth.binned <- round_any(em$Depth.m, 0.5)
  emm <- group_by(em, variable, cast, down.up, Depth.binned) %>% summarise_each(funs="mean")
  
  # pass the correct arguments to interp.dist
  my.interp.dist <- function(x, ...)  {
    xi <- interp.dist(x=x$Distance.km, y=x$Depth.m, z=x$value, ...)
    xi <- rename(xi, Distance.km=x, Depth.m=y)
    return(xi)
  }


  # interpolation across the transect
  ei <- ddply(emm, ~variable, my.interp.dist, duplicate="mean", x.step=250, y.step=1, anisotropy=1200, smooth=FALSE, theta=0.5, .progress="text")
  # TODO try kriging
  
  # compute anomalies
  eiAnom <- ddply(eiCoarse, ~variable+Depth.m, function(x) {
    mean <- mean(x$value, na.rm=TRUE)
    x$value <- x$value - mean
    return(x)
  })
  levels(eiAnom$variable) <- str_c(levels(eiAnom$variable), ".anomaly")
  eiCoarse <- rbind(eiCoarse, eiAnom)

    ggplot(mapping=aes(x=Distance.nm, y=-Depth.m, fill=value)) +
  plots <- dlply(ei, ~variable, function(Xi) {
      geom_raster(data=Xi) +
      scale_fill_spectral() + labs(title=Xi$variable[1])
  })
  pdf(str_c(dir, "/isiis_interp.pdf"), width=10, height=16)
  plots$ncol <- 2
  do.call(grid.arrange, plots)
  dev.off()
  write.csv(dcast(ei, Distance.km+Depth.m~variable, value.var="value"), file=str_c(dir, "/isiis_interp.csv"),  row.names=FALSE)

})

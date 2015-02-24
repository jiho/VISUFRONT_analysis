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
    distance <- "distanceFromStart"
  } else {
    distance <- "distanceFromShore"
  }
  
  # remove incorrect data: above 30m in downcasts
  eC <- e[-which(e$down.up %in% "down" & e$Depth.m <= 30), ]
  # ggplot(e, aes(x=distanceFromVlfr, y=-Depth.m, colour=down.up)) + geom_point()
  # ggplot(eC, aes(x=distanceFromVlfr, y=-Depth.m, colour=Salinity.PPT))  + geom_point() + scale_colour_spectral()

  # melt data to interpolate each variable sequentially
  eCm <- melt(eC, id.vars=c(distance, "Depth.m"), measure.vars=c("Salinity.PPT", "Temp.C", "Fluoro.volts", "Oxygen.ml.l", "Irrandiance.UE.cm", "Density"))
  # interp does not like NAs
  eCm <- na.omit(eCm)
  # homogenise distance name
  names(eCm)[1] <- "Distance.nm"
  
  
  # pass the correct arguments to interp.dist
  my.interp.dist <- function(x, ...)  {
    xi <- interp.dist(x=x$Distance.nm, y=x$Depth.m, z=x$value, ...)
    xi <- rename(xi, c("x"="Distance.nm", "y"="Depth.m"))
    return(xi)
  }
  
  # interpolation over a coarse grid
  eiCoarse <- ddply(eCm, ~variable, my.interp.dist, duplicate="mean", x.step=1000, y.step=2, anisotropy=1200, smooth=TRUE, theta=0.5, .progress="text")
  # TODO try kriging
  
  # compute anomalies
  eiAnom <- ddply(eiCoarse, ~variable+Depth.m, function(x) {
    mean <- mean(x$value, na.rm=TRUE)
    x$value <- x$value - mean
    return(x)
  })
  levels(eiAnom$variable) <- str_c(levels(eiAnom$variable), ".anomaly")
  eiCoarse <- rbind(eiCoarse, eiAnom)

  plots <- dlply(eiCoarse, ~variable, function(Xi) {
    ggplot(mapping=aes(x=Distance.nm, y=-Depth.m, fill=value)) +
      geom_raster(data=Xi) +
      scale_fill_spectral() + labs(title=Xi$variable[1])
  })
  pdf(str_c(dir, "/isiis_interp_coarse.pdf"), width=10, height=16)
  plots$ncol <- 2
  do.call(grid.arrange, plots)
  dev.off()
  write.csv(dcast(eiCoarse, Distance.nm+Depth.m~variable, value.var="value"), file=str_c(dir, "/isiis_interp_coarse.csv"),  row.names=FALSE)


  # # interpolation over a fine grid
  # eifine <- ddply(eCm, ~variable, my.interp.dist, duplicate="mean", x.step=250, y.step=0.5, anisotropy=1200, .progress="text")
  # # TODO try kriging
  # 
  # # compute anomalies
  # eiAnom <- ddply(eifine, ~variable+Depth.m, function(x) {
  #   mean <- mean(x$value, na.rm=TRUE)
  #   x$value <- x$value - mean
  #   return(x)
  # })
  # levels(eiAnom$variable) <- str_c(levels(eiAnom$variable), ".anomaly")
  # eifine <- rbind(eifine, eiAnom)
  # 
  # plots <- dlply(eifine, ~variable, function(Xi) {
  #   ggplot(mapping=aes(x=Distance.nm, y=-Depth.m, fill=value)) +
  #     geom_raster(data=Xi) +
  #     scale_fill_spectral() + labs(title=Xi$variable[1])
  # })
  # pdf(str_c(dir, "/isiis_interp_fine.pdf"), width=10, height=16)
  # plots$ncol <- 2
  # do.call(grid.arrange, plots)
  # dev.off()
  # write.csv(dcast(eifine, Distance.nm+Depth.m~variable, value.var="value"), file=str_c(dir, "/isiis_interp_fine.csv"),  row.names=FALSE)
  # 
})

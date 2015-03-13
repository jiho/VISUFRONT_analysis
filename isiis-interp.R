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

source("lib_plot.R")
source("lib_process.R")
data_dir <- data_dir_path()

library("reshape2")
library("plyr")
library("dplyr")
library("stringr")
library("ggplot2")
library("grid")       # for plots on a grid
library("gridExtra")
# library("doParallel")
# registerDoParallel(cores=4)

# get all isiis data
isiisFiles <- list.files(str_c(data_dir, "/transects"), pattern="isiis.csv", full=TRUE, recursive=TRUE)

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
  ei <- ddply(emm, ~variable, my.interp.dist, duplicate="mean", x.step=250, y.step=1, anisotropy=2000, smooth=TRUE, theta=0.22, .progress="text")
  # TODO try kriging
  
  # compute anomalies
  ei_anom <- group_by(ei, variable, Depth.m) %>% mutate(value=(value - mean(value, na.rm=T)))
  levels(ei_anom$variable) <- str_c(levels(ei_anom$variable), ".anomaly")
  ei <- rbind(ei, ei_anom)

  plots <- dlply(ei, ~variable, function(Xi) {
    ggplot(data=Xi, mapping=aes(x=Distance.km, y=-Depth.m, fill=value)) +
      geom_raster() +
      geom_contour(aes(z=value), colour="white", alpha=0.8, size=0.3) +
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
      # scale_fill_gradient(low="black", high="white") +
      scale_fill_spectral() +
      labs(title=Xi$variable[1]) + theme_gray(10)
  })
  pdf(str_c(dir, "/isiis_interp.pdf"), width=15, height=15)
  plots$ncol <- 2
  do.call(grid.arrange, plots)
  dev.off()
  write.csv(dcast(ei, Distance.km+Depth.m~variable, value.var="value"), file=str_c(dir, "/isiis_interp.csv"),  row.names=FALSE)
  # TODO the oxygen looks marked by the original data resolution. probably a problem between up and down casts, because of laggy sensor. Look into this.

}, .parallel=FALSE)

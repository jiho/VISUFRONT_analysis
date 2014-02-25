#
#      Useful data representation functions
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

# Spectral colour map from ColorBrewer
spectral <- function(n=6) {
  library("RColorBrewer")
  rev(brewer.pal(name="Spectral", n=n))
}

scale_fill_spectral <- function(...) {
  scale_fill_gradientn(colours=spectral(...))
}
scale_colour_spectral <- function(...) {
  scale_colour_gradientn(colours=spectral(...))
}
scale_color_spectral <- scale_colour_spectral

# Interpolate a slice of data for which the x-axis is time
interp.time <- function(x, y, z, nx=length(x)/100, y.step=2.5, ...) {
  library("akima")
  library("reshape2")

  # interpolate as numbers
  time <- as.numeric(time)
  i <- interp(x=time, y=y, z=z, xo=seq(min(time), max(time), length=nx), yo=seq(min(y), max(y), by=y.step), ...)

  # extract a data.frame
  out <- melt(i$z, varnames=c("time","y"))
  out$time <- i$x[out$time]
  out$y <- i$y[out$y]

  # reconvert to time
  out$time <- as.POSIXct(out$time, origin=as.POSIXct("1970-01-01 01:00:00.00"))

  return(out)
}

# Interpolate a slice of data for which the x-axis is a distance in nautical miles
interp.dist <- function(x, y, z, anisotropy=1000, x.step=500, y.step=2.5, smooth=FALSE, theta=0.2, ...) {
  #
  # Interpolate data over a distance coordinate
  #
  # x   vector of distance *IN NAUTICAL MILES*
  # y   vector of depth in m
  # z   vector of measured variable
  # anisotropy  anisotropy ratio between x and y
  # x/y.step    interpolation grid steps in m
  # smooth      boolean, wether to smooth the first interpolation using fields::image.smooth
  # x/y.step.smooth   interpolation grid step for the smoothing
  # grid.smooth intepolation grid for the smoothing, overrides x/y.step.smooth
  # theta       bandwidth for the kernel smoother in fields::image.smooth

  library("akima")
  library("reshape2")

  # correct x-axis for anisotropy between horizontal and vertical
  x <- x * 1852 / anisotropy

  # interpolate
  i <- interp(x=x, y=y, z=z, xo=seq(0, max(x), by=x.step/anisotropy), yo=seq(0, max(y), by=y.step), ...)

  # smooth
  if ( smooth ) {
    library("fields")
    i <- image.smooth(i, grid=list(x=i$x, y=i$y), theta=theta)
  }

  # extract a data.frame
  out <- melt(i$z, varnames=c("x","y"))
  out$x <- i$x[out$x] * anisotropy / 1852
  out$y <- i$y[out$y]

  return(out)
}

spectral <- function(n=6) {
  library("RColorBrewer")
  rev(brewer.pal(name="Spectral", n=n))
}

interp.slice <- function(x, y, z, nx=length(x)/30, ny=max(y), ...) {
  library("akima")
  library("reshape2")

  timePOSIX <- time
  time <- as.numeric(time)
  i <- interp(x=time, y=y, z=z, xo=seq(min(time), max(time), length=nx), yo=seq(min(y), max(y), length=ny), ...)
	out <- melt(i$z, varnames=c("time","y"))
	out$time <- i$x[out$time]
	out$y <- i$y[out$y]
	out <- rename(out, c(value="z"))
  
  out$time <- as.POSIXct(out$time, origin=as.POSIXct("1970-01-01 01:00:00.00"))
  
  out <- na.omit(out)
  
  return(out)
}

interp.dist <- function(x, y, z, anisotropy=1000, x.step=500, y.step=2.5, ...) {
  library("akima")
  library("reshape2")

  x <- x * 1852 / anisotropy 
  i <- interp(x=x, y=y, z=z, xo=seq(0, max(x), by=x.step/anisotropy), yo=seq(0, max(y), by=y.step), ...)
  
	out <- melt(i$z, varnames=c("x","y"))
	out$x <- i$x[out$x] * anisotropy / 1852
	out$y <- i$y[out$y]
  
  return(out)
}

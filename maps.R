#
# Produce various geographical maps to illustrate presentations
#
# (c) Copyright 2015 Jean-Olivier Irisson, GNU General Public License v3
#
#--------------------------------------------------------------------------

library("ggplot2")
library("dplyr")

# coastline, at various resolutions
coast <- read.csv("map/gshhg_coteazur.csv")
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast, fill="white")

coasth <- read.csv("map/gshhg_coteazur_h.csv")
gcoasth <- geom_polygon(aes(x=lon, y=lat), data=coasth, fill="white")

# mapping setup
setup <- list(coord_map(), scale_x_continuous(expand=c(0,0)), scale_y_continuous(expand=c(0,0)))

##{ Map of surface transects -----------------------------------------------

transects <- read.csv("ts_in_transect.csv")

transects <- filter(transects, name != "test")
range(transects$lat)
range(transects$lon)

ggplot() + geom_path(aes(x=lon, y=lat, group=name), data=transects, alpha=0.7) + gcoasth + setup

# }





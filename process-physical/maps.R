#
# Produce various geographical maps to illustrate presentations
#
# (c) Copyright 2015 Jean-Olivier Irisson, GNU General Public License v3
#
#--------------------------------------------------------------------------

source("lib_process.R")
data_dir <- data_dir_path()

library("stringr")
library("ggplot2")
library("plyr")
library("dplyr")

# coastline, at various resolutions
coast <- read.csv(str_c(data_dir, "/map/gshhg_coteazur.csv"))
gcoast <- geom_polygon(aes(x=lon, y=lat), data=coast, fill="white")

coasth <- read.csv(str_c(data_dir, "/map/gshhg_coteazur_h.csv"))
gcoasth <- geom_polygon(aes(x=lon, y=lat), data=coasth, fill="white")

# mapping setup
setup <- list(coord_map(), scale_x_continuous(expand=c(0,0)), scale_y_continuous(expand=c(0,0)))

##{ Map of surface transects -----------------------------------------------

transects <- read.csv(str_c(data_dir, "/ts_in_transects.csv"))

# remove test transect
transects <- filter(transects, name != "test")

# decompose transect name into transect type and number
transects$type <- str_replace_all(transects$name, "_[0-9]+", "")
transects$number <- str_replace_all(transects$name, "(((along_|cross_)current)|(lagrangian)|(monaco))_", "")

# make a data.frame with the starting points of each transect
starts <- ddply(transects, ~name, head, 1)

# plots
ggplot() + geom_path(aes(x=lon, y=lat, group=name), data=transects, alpha=0.7) + gcoasth + setup

ggplot() +
  geom_path(aes(x=lon, y=lat, group=name), data=transects, alpha=0.7) +
  geom_text(aes(x=lon, y=lat, label=number), data=starts, alpha=0.7, size=3, hjust=1) +
  gcoasth + setup +
  facet_wrap(~type)

# }





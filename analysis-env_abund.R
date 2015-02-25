#
#      Investigate the abundance-environment relationships
#
#  (c) Copyright 2014 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

data <- "~/Dropbox/visufront-data/"

library("stringr")
library("lubridate")
library("plyr")
library("dplyr")
library("ggplot2")
library("reshape2")
library("tidyr")
library("grid")
library("gridExtra")
library("doParallel")

source("lib_plot.R")
source("lib_zooprocess.R")

registerDoParallel(cores=4)


##{ Read abundance data from zooprocess -----------------------------------

# read identifications from zooprocess
files <- list.files(str_c(data, "/zooprocess/"), pattern=glob2rx("*_dat1.txt"), full=TRUE)
system.time(pids <- ldply(files, read.pid, .progress="text", .parallel=TRUE))
# clean names
(taxa <- sort(unique(pids$Valid)))

# remove unwanted, numerous particles
d <- filter(pids, !(Valid %in% c("noise", "bad_focus", "det_aggregates", "det_fibers", "det_aggregates_day", "duplicates", "unidentified", "unidentified_of_interest")))

# basic counts
# number of organisms
nrow(d)
# [1] 113044
# number of groups
length(unique(d$Valid))

# TODO group some taxa together

# compute abundance per depth bin
bin <- 1
d$depth <- round_any(d$Depth, bin)
d <- group_by(d, Valid, Label, depth) %>% summarise(Abund=n())

# add profile number
label_bits <- str_split_fixed(d$Label, fixed("_"), 3)
d$profile <- as.numeric(label_bits[, 3])
# add transect number
d$transect <- label_bits[, 2]

# check abundances
print(d %>% group_by(transect, Valid) %>%
    summarise(tot_abund=sum(Abund)) %>%
    spread(key=transect, value=tot_abund), n=50)
# TODO one NA group in cc4. check this
d <- filter(d, !is.na(Valid))

# get volume sampled per bin
# get frame (2048x2048 image) names from the datfile
files <- list.files(str_c(data, "/zooprocess/"), pattern=glob2rx("*_datfile.txt"), full=TRUE)
dat <- adply(files, 1, read.table, sep=";", strip.white=TRUE, .progress="text")
# bin the depth over the same bins
dat$depth <- round_any(dat$V3/10, bin)
# add transect+profile Label
file_names <- basename(files)
labels <- str_replace(file_names, "_datfile.txt", "")
dat$Label <- labels[dat$X1]
# compute volume from number of frames
dat <- group_by(dat, Label, depth) %>% summarise(n_frames=n())
dat$vol.m3 <- dat$n_frames * 7.78 / 1000


# get volume in the biological data
d <- left_join(d, select(dat, -n_frames))
# and compute concentration
d$abund.m3 <- d$Abund / d$vol.m3

# cleanup data.frame
d <- select(ungroup(d), transect, profile, depth, taxon=Valid, abund.m3)

# }
# save(d, bin, file="isiis_catches.Rdata")


##{ Geolocalize abundance data ---------------------------------------------

# load("isiis_catches.Rdata")
# get distance measure, depth, and date/time from the hydrological data
# this will allow us to extract env data from the interpolated record

# read raw hydrological data
h4 <- read.csv("transects/cross_current_4/isiis.csv")
h4$transect <- "cc4"
h5 <- read.csv("transects/cross_current_5/isiis.csv")
h5$transect <- "cc5"
h <- rbind(h4, h5)

# keep only down casts (the biological data is on downcasts only)
h <- filter(h, down.up=="down")

# keep only the transects x cast investigated in the biological data
h <- filter(h, interaction(transect, profile) %in% unique(interaction(d$transect, d$profile)))

# parse date and time
options(digits.secs=2)
h$dateTimeMsec <- ymd_hms(h$dateTimeMsec)
h$dateTime <- ymd_hms(h$dateTime)

# bin depth over the same bins as biological data to be able to join them
h$depth <- round_any(h$Depth.m, bin)

# average location data per depth bin
loc <- h %>%
      select(transect, profile, depth, dateTime=dateTimeMsec, starts_with("dist"), lat, lon) %>%
      group_by(transect, profile, depth) %>%
      summarise_each(fun="mean")

# some depth bins are not covered; interpolate the location data in those bins
loc <- ddply(loc, ~transect+profile, function(x) {

  obs_depths <- x$depth                                  # observed
  all_depths <- seq(min(x$depth), max(x$depth), by=bin)  # full range

  # linear interpolation
  linint <- function(y, x_from=obs_depths, x_to=all_depths) {
    approx(x=x_from, y=y, xo=x_to)$y
  }
  
  # interpolate all data
  dateTime        <- as.POSIXct(linint(x$dateTime), origin=ymd("1970-01-01"), tz="UTC")
  dist_from_start <- linint(x$dist_from_start)
  dist_from_vlfr  <- linint(x$dist_from_vlfr)
  dist_from_shore <- linint(x$dist_from_shore)
  lat             <- linint(x$lat)
  lon             <- linint(x$lon)
  
  return(data.frame(depth=all_depths, dateTime, dist_from_start, dist_from_vlfr, dist_from_shore, lat, lon))
})

# these are all the actually sampled locations/times
# in the biological data, we only have data where something was recorded
# joint the two to add zeros in the biological record when nothing was captured
# compute all possibilities of location x taxa
taxa <- sort(unique(d$taxon))
all_loc <- ddply(loc, .(transect, profile, depth), function(x) {
  data.frame(x, taxon=taxa)
})

# for all points sampled, add captures where we have some
all <- left_join(all_loc, d)
# this puts NAs where there are no captures
# replace NAs by 0
all$abund.m3[is.na(all$abund.m3)] <- 0

# check how many geolocalised capture we cannot match with hydrological data
# there should not be any... but there are because of differences in rounding depth probably
allf <- full_join(all_loc, d)
nrow(allf) - nrow(all)
group_by(allf, transect, profile) %>% summarise(sum(is.na(lat)))
unique(filter(allf, is.na(lat))[,1:3])
# -> OK, just a few missing bits, forget about it

# convert to wide format for multivariate stats
d <- all
dw <- spread(d, key=taxon, value=abund.m3)

# }
save(d, dw, file="isiis_catches.Rdata")


##{ Associate environmental data with biological records ------------------

load("isiis_catches.Rdata")
# NB: we read environmental data from the interpolated version to
#     - smooth out small scale variability
#       (that might not be a good idea for fine scale patterns however...)
#     - use anomalies for env variables, which are only computed on the interpolated fields

# read interpolated data
ei4 <- read.csv("transects/cross_current_4/isiis_interp.csv")
ei4$transect <- "cc4"
ei5 <- read.csv("transects/cross_current_5/isiis_interp.csv")
ei5$transect <- "cc5"
ei <- rbind(ei4, ei5)


# Extract environmental data in X to locations in xy
get.env <- function(env, locs) {
  # prepare a x,y,z list (persp-like) for interp.surface
  m <- acast(env, Distance.km ~ Depth.m, value="value")
  x <- sort(unique(env$Distance.km))
  y <- sort(unique(env$Depth.m))
  xyz <- list(x=x, y=y, z=m)
  # image(xyz)
  library("fields")
  i <- interp.surface(xyz, loc=as.matrix(locs))
  return(data.frame(locs, value=i))
}

env_at_loc <- ddply(dw, ~transect, function(x) {
  # get locations of biological data for this transect
  this_loc <- select(x, dist_from_shore, depth)

  # get interpolated environmental data for this transect
  this_env <- filter(ei, transect == x$transect[1])
  # convert it to tall format to treat each variable separately
  this_env <- gather(this_env, key=variable, value=value, -transect, -Distance.km, -Depth.m)
  
  # extract env data (for each variable) at those locations
  env_at_loc <- group_by(this_env, transect, variable) %>% do(get.env(., locs=this_loc))

  # convert back to wide format
  env_at_loc <- spread(env_at_loc, key=variable, value=value)

  return(env_at_loc)
})

# and environmental data to the geolocalised captures
d <- left_join(d, env_at_loc)
dw <- left_join(dw, env_at_loc)

# }
save(d, dw, file="isiis_catches.Rdata")


##{ Regression trees ------------------------------------------------------

plot.abund.env.map <- function(taxon, variable) {
  ggplot() + geom_raster(aes_string(x="Distance.nm", y="-Depth.m", fill=variable), data=ei) + scale_fill_spectral() + geom_point(aes_string(x="distanceFromVlfr", y="-DepthBin", size=taxon), data=D) + scale_size_area() + scale_x_continuous(limits=c(0,30))
}

plot.abund.env <- function(taxon, variable) {  
  d <- D[which(D[,taxon] > 0),]
  d[,taxon] <- log1p(d[,taxon])
  ggplot(d, aes_string(x=variable, y=taxon)) + geom_point() + geom_smooth()
}

fit.tree <- function(taxon, ...) {
  library("stringr")
  d <- D[which(D[,taxon] > 0),]
  d[,taxon] <- log1p(d[,taxon])
  formula <- str_c(taxon , "~ Fluoro.volts + Oxygen.ml.l + Oxygen.ml.l.anomaly + Salinity.PPT + Salinity.PPT.anomaly + Temp.C.anomaly")
  library("mvpart")
  mvpart(formula, data=d, ...)
}

m <- fit.tree("doliolids", xv="1se", xvmult=100, xval=5)
plot.abund.env.map("doliolids", "Oxygen.ml.l")
plot.abund.env.map("doliolids", "Fluoro.volts")

m <- fit.tree("ephyrae", xv="1se", xvmult=100, xval=5)
plot.abund.env.map("ephyrae", "Oxygen.ml.l")
plot.abund.env.map("ephyrae", "Salinity.PPT")
plot.abund.env("ephyrae", "Oxygen.ml.l")
plot.abund.env("ephyrae", "Salinity.PPT")

m <- fit.tree("fish", xv="1se", xvmult=100, xval=5)
plot.abund.env.map("fish", "Temp.C.anomaly")
plot.abund.env.map("fish", "Salinity.PPT")

m <- fit.tree("radiolarian_dark", xv="1se", xvmult=100, xval=5)
plot.abund.env.map("radiolarian_dark", "Fluoro.volts")
plot.abund.env.map("radiolarian_dark", "Oxygen.ml.l")
plot.abund.env("radiolarian_dark", "Fluoro.volts")

# }


##{ Boosted regression trees ----------------------------------------------

# library("gbm")
# 
# fit.tree <- function(taxon, verbose=FALSE, ...) {
#   library("stringr")
#   D <- dd[which(dd[,taxon] > 0),]
#   D[,taxon] <- round(D[,taxon] * 100)
#   formula <- as.formula(str_c(taxon , " ~ Fluoro.volts + Oxygen.ml.l + Oxygen.ml.l.anomaly + Salinity.PPT + Salinity.PPT.anomaly + Temp.C.anomaly"))
#   m <- gbm(
#     formula=radiolarian_dark ~ Fluoro.volts + Oxygen.ml.l + Oxygen.ml.l.anomaly + Salinity.PPT + Salinity.PPT.anomaly + Temp.C.anomaly,
#     data=D,
#     distribution="poisson",
#     n.trees=1000,         # should be > 1000 to be robust
#     shrinkage=0.001,      # should be small to allow enough trees
#     cv.folds=5,           # allows estimating error in prediction and then use gbm.perf with method cv to find optimal n.trees
#     interaction.depth=4,  # higher level interactions means faster fit => decrease shrink/learning rate to compensate
#     bag.fraction=0.5,     # proportion of data randomly drawn to compute each tree, adds robustness to the model. Allows to use gbm.perf with method OOB to find optimal n.trees
#     verbose=FALSE
#   )
# }
# 
# m <- fit.tree(taxon="radiolarian_dark")
#
# }

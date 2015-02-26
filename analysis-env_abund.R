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
save(d, bin, file="zooprocess_data.Rdata")


load("zooprocess_data.Rdata")
##{ Geolocalize abundance data ---------------------------------------------

# get distance measure, depth, and date/time from the hydrological data
# this will allow us to extract env data from the interpolated record

# read raw hydrological data
h4 <- read.csv("transects/cross_current_4/isiis.csv")
h4$transect <- "cc4"
h5 <- read.csv("transects/cross_current_5/isiis.csv")
h5$transect <- "cc5"
h <- rbind(h4, h5)

# total number of casts
nrow(unique(select(h, transect, cast))

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
# there should not be any...
allf <- full_join(all_loc, d)
nrow(allf) - nrow(all)
# -> OK!
# group_by(allf, transect, profile) %>% summarise(sum(is.na(lat)))
# unique(filter(allf, is.na(lat))[,1:3])

# convert to wide format for multivariate stats
d <- all
dw <- spread(d, key=taxon, value=abund.m3)

# }
save(d, dw, file="zooprocess_data_geoloc.Rdata")


load("zooprocess_data_geoloc.Rdata")
##{ Distribution of a few groups -------------------------------------------

# flag day and night
dw$time_of_day <- factor(dw$transect, levels=c("cc4", "cc5"), labels=c("night", "day"))
d$time_of_day <- factor(d$transect, levels=c("cc4", "cc5"), labels=c("night", "day"))

# read interpolated data
ei4 <- read.csv("transects/cross_current_4/isiis_interp.csv")
ei4$transect <- "cc4"
ei5 <- read.csv("transects/cross_current_5/isiis_interp.csv")
ei5$transect <- "cc5"
ei <- rbind(ei4, ei5)
ei$time_of_day <- factor(ei$transect, levels=c("cc4", "cc5"), labels=c("night", "day"))

# plot physics
base <- ggplot(data=ei, aes(x=Distance.km, y=-Depth.m)) +
          facet_grid(time_of_day~.) +
          labs(x="Distance from shore (km)", y="Depth (m)") +
          scale_fill_spectral(na.value=NA) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
sal_contour <- geom_contour(aes(z=Salinity.PPT), breaks=38.25, colour="white", na.rm=T, alpha=0.7)
base + geom_raster(aes(fill=Salinity.PPT))
base + geom_raster(aes(fill=Salinity.PPT)) + sal_contour
base + geom_raster(aes(fill=Fluoro.volts)) + sal_contour
base + geom_raster(aes(fill=Temp.C)) + sal_contour
base + geom_raster(aes(fill=Temp.C.anomaly)) + sal_contour + labs(fill="Temp.C\nanomaly")
base + geom_raster(aes(fill=Density.anomaly)) + sal_contour

# pot biological groups
# base plot
base <- ggplot(data=dw, aes(x=dist_from_shore, y=-depth)) +
          geom_point(size=0.4, alpha=0.2) +
          facet_grid(time_of_day~.) +
          geom_contour(aes(x=Distance.km, y=-Depth.m, z=Salinity.PPT), data=ei, breaks=38.25, colour="grey50", na.rm=T) +
          scale_size_area() +
          labs(x="Distance from shore (km)", y="Depth (m)")

(taxa <- unique(d$taxon))

base + geom_point(aes(size=fish_larvae)) 
# surface, coastal
base + geom_point(aes(size=doliolids))
# surface, coastal
base + geom_point(aes(size=phyto_diatom_chains))
# dcm?, offshore
base + geom_point(aes(size=append_fritill))
# deep
base + geom_point(aes(size=chaetognaths))
# ~surface, coastal
base + geom_point(aes(size=cteno_mertens))
# coastal, deep
base + geom_point(aes(size=ephyrae))
# coastal strong migration
base + geom_point(aes(size=radiolarians_solitarian_dark)) + labs(size="radiolarian\nsolitary\ndark")


# }
save(d, dw, ei, file="zooprocess_data_env.Rdata")

load("zooprocess_data_env.Rdata")
##{ Associate environmental data with biological records ------------------

# NB: we read environmental data from the interpolated version to
#     - smooth out small scale variability
#       (that might not be a good idea for fine scale patterns however...)
#     - use anomalies for env variables, which are only computed on the interpolated fields

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
  this_env <- gather(this_env, key=variable, value=value, -transect, -Distance.km, -Depth.m, -time_of_day)
  
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
save(d, dw, file="zooprocess_data_with_env.Rdata")

load("zooprocess_data_with_env.Rdata")
##{ Regression trees ------------------------------------------------------

fit.tree <- function(taxon, ...) {
  library("stringr")
  # select taxon of interest
  this_d <- dw[which(dw[,taxon] > 0),]
  # log-transform catches to smooth variance
  this_d[,taxon] <- log1p(this_d[,taxon])
  # create model
  formula <- str_c(taxon , "~ Fluoro.volts + Fluoro.volts.anomaly + Oxygen.ml.l + Oxygen.ml.l.anomaly + Salinity.PPT + Salinity.PPT.anomaly + Temp.C + Temp.C.anomaly")
  # fit and inspect model
  library("mvpart")
  m <- mvpart(formula, data=this_d, xv="1se", xvmult=100, xval=5, ...)
  # NB: mvpart is deprecated and we don't use the multivariate aspect anyway so we should reimplement it with rpart
  # library("rpart")
  # m <- rpart(formula, data=this_d, rpart.control(minbucket=5, xval=5, cp=0.02))
  # plot(m)
  # text(m)
  # print(m$variable.importance)

  return(m)
}

plot.abund.env.map <- function(taxon, variable) {
  ggplot(data=dw, aes(x=dist_from_shore, y=-depth)) +
    geom_raster(aes_string(x="Distance.km", y="-Depth.m", fill=variable), data=ei) +
    geom_contour(aes(x=Distance.km, y=-Depth.m, z=Salinity.PPT), data=ei, breaks=38.25, colour="white", na.rm=T, alpha=0.7) +
    facet_grid(time_of_day~.) +
    geom_point(size=0.4, alpha=0.2) +
    scale_fill_spectral(na.value=NA) +
    geom_point(aes_string(size=taxon)) +
    scale_size_area() +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    labs(x="Distance from shore (km)", y="Depth (m)")
}

(taxa <- unique(d$taxon))

m <- fit.tree("doliolids")
# -> temp + salinity
plot.abund.env.map("doliolids", "Temp.C")
plot.abund.env.map("doliolids", "Salinity.PPT")

m <- fit.tree("ephyrae")
# bof

m <- fit.tree("fish_larvae")
# bof

m <- fit.tree("radiolarians_solitarian_dark")
# -> fluo, CV error 0.356
plot.abund.env.map("radiolarians_solitarian_dark", "Fluoro.volts") + labs(size="radiolarian\nsolitary\ndark")

# }


##{ Boosted regression trees ----------------------------------------------

fit.brt <- function(taxon, verbose=FALSE, ...) {
  library("stringr")
  dd <- dw[which(dw[,taxon] > 0),]
  dd[,taxon] <- round(dd[,taxon] * 10)
  vars <- c("Fluoro.volts", "Fluoro.volts.anomaly", "Oxygen.ml.l", "Oxygen.ml.l.anomaly", "Salinity.PPT", "Salinity.PPT.anomaly", "Temp.C", "Temp.C.anomaly")
  dd <- dd[,c(taxon, vars)]
  formula <- as.formula(str_c(taxon , " ~ ."))
  library("gbm")
  m <- gbm(
    formula=formula,
    data=dd,
    distribution="poisson",
    ...
  )
  print(m)
  print(summary(m))
  return(m)
}

# fit model
m <- fit.brt(taxon="radiolarians_solitarian_dark", n.trees=4000, shrinkage=0.02, cv.folds=5, interaction.depth=4, verbose=FALSE)

# predict distribution from environmental variables
ei_var <- select(ei, -Distance.km, -Depth.m, -transect, -time_of_day)
p <- predict(m, newdata=ei_var, n.trees=gbm.perf(m, method="cv", plot.it=FALSE), type="response")
er <- cbind(ei, abund.m3=p/10)
ggplot(data=er, aes(x=Distance.km, y=-Depth.m)) +
  geom_raster(aes(fill=abund.m3)) +
  geom_contour(aes(x=Distance.km, y=-Depth.m, z=Salinity.PPT), data=ei, breaks=38.25, colour="white", na.rm=T, alpha=0.7) +
  facet_grid(time_of_day~.) +
  scale_fill_bw(na.value=NA) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  labs(x="Distance from shore (km)", y="Depth (m)")

# }

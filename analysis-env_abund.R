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
library("grid")
library("gridExtra")
library("doParallel")

source("lib_plot.R")
source("lib_zooprocess.R")

registerDoParallel(cores=4)


##{ Read abundance data from zooprocess -----------------------------------

# read identifications from zooprocess
files <- list.files(str_c(data, "/zooprocess/"), pattern=glob2rx("*_dat1.txt"), full=TRUE)
pids <- ldply(files, read.pid, .progress="text", .parallel=TRUE)
# clean names
pids$Valid[which(pids$Valid %in% c("sipho_tail", "sipho_round"))] <- "sipho"
pids$Valid <- str_replace(pids$Valid, "ephyrae_side", "ephyrae")
pids$Valid <- str_replace(pids$Valid, "radiolarians", "radiolarian")
pids$Valid <- str_replace(pids$Valid, "shrimp_large", "shrimps")
pids$Valid <- str_replace(pids$Valid, "polychaets", "polychaetes")

# TODO group some taxa together

# compute abundance per depth bin
bin <- 1
pids$Depth_binned <- round_any(pids$Depth, bin)
bio <- group_by(pids, Valid, Label, Depth_binned) %>% summarise(Abund=n())

# add cast number
label_bits <- str_split_fixed(bio_all$Label, fixed("_"), 3)
bio_all$cast <- as.numeric(label_bits[, 3]) * 2 - 1
# ggplot(bio_all[bio_all$Valid == "fish",]) + geom_point(aes(y=-Depth_binned, x=Label, size=Abund)) # + scale_size_area()
# add transect number
bio_all$transect <- label_bits[, 2]

# get volume sampled per bin
# count number of frames from the datfile
files <- list.files(str_c(data, "/zooprocess/"), pattern=glob2rx("*_datfile.txt"), full=TRUE)
dat <- adply(files, 1, read.table, sep=";", strip.white=TRUE, .progress="text")
dat$Depth_binned <- round_any(dat$V3/10, bin)

# add transect name
fileNames <- basename(files)
transects <- str_replace(fileNames, "_datfile.txt", "")
dat$Label <- transects[dat$X1]

# compute volume from number of frames
vol <- group_by(dat, Label, Depth_binned) %>% summarise(n_frames=n())
vol$vol.m3 <- vol$n_frames * 7.78 / 1000


# get volume in the biological data
d <- left_join(bio_all, select(vol, -n_frames))
# and compute concentration
d$abund.m3 <- d$Abund / d$vol.m3

# cleanup data.frame
d <- select(d, transect, cast, depth=Depth_binned, taxon=Valid, abund=Abund, abund.m3)

# convert to wide format for multivariate stats
dd <- dcast(d, transect+cast+depth~taxon, value.var="abund.m3")

# }
save(d, dd, file="isiis_catches.Rdata")

##{ Associate environmental data with biological records ------------------

# read env data
e <- read.csv("transects/cross_current_4/isiis.csv")
ei <- read.csv("transects/cross_current_4/isiis_interp_coarse.csv")

# get distance associated with biological data to be able to extract environmental data associated with biological data
# compute distances for each depth bin
e$DepthBin <- round_any(e$Depth.m, bin)
eLoc <- ddply(e, ~cast+down.up+DepthBin, function(x) {
  colMeans(x[,c("lat", "lon", "distanceFromStart", "distanceFromVlfr", "distanceFromShore")])
}, .progress="text")
# get the distances
dLoc <- join(dd, eLoc, type="inner")

# interpolate env data at locations of biological data
xy <- dLoc[,c("distanceFromVlfr", "DepthBin")];

get.env <- function(X, xy) {
  m <- acast(X, Distance.nm ~ Depth.m, value="value")
  x <- sort(unique(X$Distance.nm))
  y <- sort(unique(X$Depth.m))
  xyz <- list(x=x, y=y, z=m)
  # image(xyz)
  library("fields")
  i <- interp.surface(xyz, loc=xy)
  return(i)
}

eim <- melt(ei, id.vars=c("Distance.nm", "Depth.m"))
eiB <- ddply(eim, ~variable, get.env, xy=xy)
eiBd <- as.data.frame(t(eiB[,-1]))
names(eiBd) <- eiB[,1]

D <- cbind(dLoc, eiBd)

# }


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

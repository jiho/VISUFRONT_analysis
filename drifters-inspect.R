# Check relative positions and time

source("lib_process.R")
data_dir <- data_dir_path()

library("gdata")
library("lubridate")
library("stringr")
library("plyr")
library("dplyr")
library("ggplot2")


# read drifter trajectories
d <- read.xls(str_c(data_dir, "/drifters/drifters.xls"), na.strings="", stringsAsFactors=FALSE)
b <- read.xls(str_c(data_dir, "/drifters/boa.xls"), na.strings="", stringsAsFactors=FALSE)
f <- read.xls(str_c(data_dir, "/drifters/float.xls"), na.strings="", stringsAsFactors=FALSE)

# compute lat and lon for boa
latBits <- str_split_fixed(b$lat, fixed("."), 2)
b$lat <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
b$lat <- 43 + b$lat/60
lonBits <- str_split_fixed(b$lon, fixed("."), 2)
b$lon <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
b$lon <- 7 + b$lon/60

# combine all data
d <- rbind(d,b,f)

# d <- filter(d, unit != "provbio")

# sort by date and time
d$dateTime <- ymd_hms(str_c(d$date, d$time, sep=" "))
d <- arrange(d, dateTime)

# plot and display time
ggplot(head(d, 80)) + geom_point(aes(x=lon, y=lat, colour=unit)) + geom_text(aes(x=lon, y=lat, colour=unit, label=dateTime), size=2, hjust=1.1)

# -> all are on the same time reference: UTC
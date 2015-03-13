#
#
# (c) Copyright 2014 Jean-Olivier Irisson, GNU General Public License v3
#
#--------------------------------------------------------------------------

source("lib_process.R")
# data_dir <- data_dir_path()

library("stringr")
library("lubridate")
library("plyr")
library("dplyr")


##{ Read all N1R files -----------------------------------------------------

# list N1R files
n1r_files <- list.files(str_c(data_dir, "/_raw_/ADCP172_000000.LTA/"), pattern=glob2rx("*.N1R"), full=TRUE)

# read all files
n1r <- ldply(n1r_files, function(f) {
  x <- scan(f, what="character", sep="\n", quiet=TRUE)
  data.frame(text=x, file=basename(f))
}, .progress="text")

# }

# dates
(dates <- ddply(n1r, ~file, function(x) {
  # detect ADCP records
  x <- x$text[str_detect(x$text, "^\\$PADCP")]
  n <- length(x)
  x <- x[str_length(x) %in% 34:40]
  n_clean <- length(x)
  # parse date and time
  x <- read.table(text=paste(x, collapse="\n"), sep=",", colClasses="character")
  date_time <- ymd_hms(str_c(x$V3, " " ,x$V4), tz="UTC")
  date_time <- with_tz(date_time, tz="Europe/Paris")

  # data.frame(n, n_clean, start=min(date_time), end=max(date_time))
  data.frame(n, start=min(date_time), end=max(date_time))
}))


# Number of correct/missing HDT
(hdt <- ddply(n1r, ~file, function(x) {
  x <- x$text[str_detect(x$text, "^\\$GPHDT") & !str_detect(x$text, "\\$PADCP")]
  n <- length(x)

  b <- read.table(text=paste(x, collapse="\n"), sep=",")
  n_missing <- sum(is.na(b$V2) | (b$V2 < 0 | b$V2 > 360))
  return(data.frame(n, n_missing, prop_ok=(1-(n_missing/n))*100))
}))

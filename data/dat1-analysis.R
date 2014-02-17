# -----------------------------------------------------------
#
# Process ISIIS data for the Ocean Science Meeting 2014
# Robin Faillettaz, 2014-02-14
#
# -----------------------------------------------------------
#
# Working on transect cross - front 04
# It's a night transect
# Profil validated : 01, 03, 05, 09, 13, 17, 25 (?)
#
# A .Rprofile file should be set personnaly to the directory where the data are stored (dropbox)

# For now locate data from the dropbox repo manually
dir <- "/Users/faillettaz/Dropbox/robin/visufront-data/"


library("plyr")
library("ggplot2")
source("~/r-utils/lib_zooprocess.R")


files <- list.files(dir, full = T)
pid <- read.pid(files[2])



pids <- NA
pids <- adply(files, 1, function(files)
    {
    pid <- read.pid(files)
    # make validation column name uniform
    names <- names(pid)
    #names <- str_replace_all(names, "Valid", "Valid")
    # force latest validation to have the name "Valid"
    validationColumns <- which(str_detect(names, "Valid"))
    if (length(validationColumns) == 1) {
         names[validationColumns] <- "Valid"
    }
    if (length(validationColumns) > 1) {
        names[max(validationColumns)] <- "Valid"    
    }
   
    names(pid) <- names
    rbind(pids, pid)
    }, .progress="text")  
pids <- pids[-1, -1]

head(pids)

# set depth BIN 
pids$DepthBin <- round_any(pids$Depth, 1)

# get raw abundance
d <- ddply(pids, ~Valid+Label+DepthBin, function(x) {sum(na.omit(x$Valid==paste(x$Valid)))})

d <- rename(d, c("V1" = "Abund"))

# plot everything
ggplot(d) + 
    geom_point(aes(x = Abund, y = -DepthBin, colour = Label)) + 
    geom_path(aes(x = Abund, y = -DepthBin, colour = Label)) + 
    facet_grid(.~Valid, scales="free_x")



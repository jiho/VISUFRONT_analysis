
library("stringr")

source("lib_zooprocess.R")  # for read.pid
source("lib_confusion.R")   # for confusion_matrix, confusion_stats, etc.

# prefix for the data location
data <- "C:/Users/petitesorciere/Dropbox/visufront-data/"

# read data
# combine info from several transects
# read data with the same categories predicted
d01 <- read.pid(str_c(data, "zooprocess/transect_cc4_01_dat1.txt"))
d03 <- read.pid(str_c(data, "zooprocess/transect_cc4_03_dat1.txt"))
d05 <- read.pid(str_c(data, "zooprocess/transect_cc4_05_dat1.txt"))
d17 <- read.pid(str_c(data, "zooprocess/transect_cc4_17_dat1.txt"))
d25 <- read.pid(str_c(data, "zooprocess/transect_cc4_25_dat1.txt"))

d <- rbind(d01,d03,d05,d17,d25)


# reduce to categories of the learning set
lg <- read.csv(str_c(data, "/zooprocess/learning_groups.csv"))
d$ValidReduced <- lg$prediction[match(d$Valid, lg$validation)]


# compute confusion matrix
cm <- confusion_matrix(d$Pred, d$ValidReduced)


# plot confusion matrix
autoplot(cm)
autoplot(cm, norm="row")
autoplot(cm, norm="col")

# compute confusion statistics (recall, precision, ...)
confusion_stats(cm)
confusion_stats(cm, sort.by="recall")
confusion_stats(cm, sort.by="precision")
confusion_stats(cm, sort.by="F1")

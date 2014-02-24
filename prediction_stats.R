
library("stringr")

source("lib_zooprocess.R")  # for read.pid
source("lib_confusion.R")   # for confusion_matrix, confusion_stats, etc.

# prefix for the data location
data <- "~/Dropbox/visufront-data/"

# read data
d <- read.pid(str_c(data, "zooprocess/transect_cc4_01_dat1.txt"))

# compute confusion matrix
cm <- confusion_matrix(d$Pred, d$Valid)

# plot confusion matrix
autoplot(cm)
autoplot(cm, norm="row")
autoplot(cm, norm="col")

# compute confusion statistics (recall, precision, ...)
confusion_stats(cm)
confusion_stats(cm, sort.by="recall")
confusion_stats(cm, sort.by="precision")
confusion_stats(cm, sort.by="F1")

# reduce to categories of the learning set
lg <- read.csv(str_c(data, "/zooprocess/learning_groups.csv")
d$ValidReduced <- lg$prediction[match(d$Valid, lg$validation)]

cm <- confusion_matrix(d$Pred, d$ValidReduced)
autoplot(cm)
autoplot(cm, "row")
autoplot(cm, "col")
confusion_stats(cm)

# combine info from several transects
# rbind()


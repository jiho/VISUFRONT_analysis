
library("stringr") # use for str
library("FactoMineR") # use for PCA
source("lib_zooprocess.R")
source("lib_confusion.R")

# read data
data <- "C:/Users/petitesorciere/Dropbox/visufront-data/"
isiis<-read.table("isiis.csv",header=TRUE,sep=",",dec=".")

# PCA, variables factor map
pcaisiis<-PCA(isiis[,2:8])
round(pcaisiis$eig,2)

# what explain axes 
barplot(pcaisiis$eig[,1])
lapply(dimdesc(pcaisiis),lapply,round,2)

dist<-dist(isiis[,2:8],method="euclidean") # memory.size problem

# relationship Depth-abundance
# read data
library("plyr")
d01 <- read.pid(str_c(data, "zooprocess/transect_cc4_01_dat1.txt"))
d03 <- read.pid(str_c(data, "zooprocess/transect_cc4_03_dat1.txt"))
d05 <- read.pid(str_c(data, "zooprocess/transect_cc4_05_dat1.txt"))
d17 <- read.pid(str_c(data, "zooprocess/transect_cc4_17_dat1.txt"))
d25 <- read.pid(str_c(data, "zooprocess/transect_cc4_25_dat1.txt"))

d <- rbind(d01,d03,d05,d17,d25)

d$DepthRound <- round_any(d$Depth, 2)
# Plot Valid=f(Depth)
ggplot(data = d)  + geom_point(aes(x= Valid , y= -Depth), alpha = 0.01)
ggplot(data = d)  + geom_boxplot(aes(x= Valid , y= -Depth))
ggplot(data = d)  + geom_violin(aes(x= Valid , y= -Depth))


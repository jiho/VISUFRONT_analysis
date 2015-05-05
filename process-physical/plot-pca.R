
plot_pca <- function(x, colour = 1){

varNames <- rownames(x$var$coord)
indNames <- rownames(x$ind$coord)

varCoord <- data.frame(x=x$var$coord[,1], y=x$var$coord[,2], size=(x$var$cos2[,1]))
indCoord <- data.frame(x=x$ind$coord[,1], y=x$ind$coord[,2], size=x$ind$cos2[,1], colour=as.factor(colour))#, shape=as.factor(shape))

VarAxe1 <- round(x$eig[1,2], digits=2)
VarAxe2 <- round(x$eig[2,2], digits=2)

varContrib <- round(x$var$contrib[,1], digits=2)
indContrib <- round(x$ind$contrib[,1], digits=2)

ggplot()+
  geom_point(aes(x=x, y=y, size=size, colour=colour), data=indCoord)+ #shape=c(rep(1, times=17), rep(17, times=17))
  #geom_text(aes(x=x, y=y,label=indNames,  color=year), size=4, data=indCoord, vjust=1.5)+
  #geom_point(aes(x=x, y=y, size=size), data=varCoord, colour="gray2", shape=17)+
  geom_text(aes(x=x, y=y, label=varNames), data=varCoord,  size=4, vjust=1.5, colour="gray2")+
  geom_segment(aes(x=0, xend=x, y=0, yend=y), arrow=arrow(length= unit(0.2,"cm"), type="closed"), data=varCoord)+
  geom_vline(aes(x=0), linetype=2, size=0.6)+  
  geom_hline(aes(y=0), linetype=2, size=0.6)+ 
  scale_x_continuous(name=paste("Axe1 (", VarAxe1, "%)", sep=" "))+
  scale_y_continuous(name=paste("Axe2 (", VarAxe2, "%)", sep=" "))+
  scale_size_continuous(paste("cos2"), guide="none") +
  scale_colour_discrete("") +
  #scale_shape_manual(name="", values=c(1, 17, 25, 4), guide = "none") + 
  opts


}
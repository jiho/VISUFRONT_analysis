#--------------------------------------------------------------#
#
#		Process and plots of the plankton data from
#		the VISUFRONT cruise in July 2013
#		Robin Faillettaz - 2013/12/10
# 
#--------------------------------------------------------------#


library("ggplot2")
library("plyr")
library("stringr")
library("reshape2")

source("data/lib_process.R")



# read coastline to plot the trajectories and check
coast <- read.csv("map/cote_azur.csv")
load("map/coast_bathy.RData")

# read stations position
station <- read.csv("Plankton-nets/station-regent-visufront.csv", header=T, sep=";")
head(station)

# compute lat and lon for stations
latBits <- str_split_fixed(station$lat_out, fixed("."), 2)
station$lat <- as.numeric(latBits[,1]) + (as.numeric(latBits[,2])/60)
station$lat <- 43 + station$lat/60
lonBits <- str_split_fixed(station$lon_out, fixed("."), 2)
station$lon <- as.numeric(lonBits[,1]) + (as.numeric(lonBits[,2])/60)
station$lon <- 7 + station$lon/60

# Add Boussole as the reference point
Boussole <- data.frame(lat=43.38, lon=7.83)

p <- ggplot(mapping=aes(x=lon, y=lat)) + 
	geom_polygon(data=coast, fill="grey60") + 
	geom_point(aes(color=date), data=station, size=3) + 
	geom_text(data=station, label=station$station_nb, vjust=2) + 
	geom_point(data=Boussole, color="blue") + 
	coord_map(xlim=c(7,8.1), ylim=c(43.2,43.75))
print(p)



# plot boat traj + stations (from drifter-plot.R)
# read ship trajectory from ts
filenames <- list.files("TS/")

s <- adply(filenames, 1, function(x) {
	s <- read.ts(str_c("TS/",x))
	return(s)
	})

ggplot(mapping=aes(x=lon, y=lat)) + 
	#geom_raster(aes(fill=-z, x=x, y=y), data=bathyDF) + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) + 
	geom_polygon(fill="gray25", data=coast, aes(x=lon, y=lat)) +
	geom_path(size=0.4, na.rm=T, data=s) + # ship track 
	geom_point(aes(color=date), data=station, size=3) + 
	geom_text(data=station, label=station$station_nb, vjust=2) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="blue") +
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="blue", vjust=2, hjust=1, size=4) +
	scale_x_continuous("Longitude", expand=c(0,0)) + 
	scale_y_continuous("Latitude", expand=c(0,0)) +
	#scale_fill_gradient(name="Depth") +
	scale_color_discrete("") + 
	coord_quickmap(xlim=c(6.8, 8.1), ylim=c(43.2, 43.75)) +
	theme_bw()



###--------  LARVAL FISH ABUNDANCE AND DISTRIBUTION  -------###


# for fish eggs 
ggplot(mapping=aes(x=lon, y=lat)) + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(data=coast, fill="grey60") +  
	geom_point(aes(color=date, size=eggs), data=station) + 
	geom_point(data=Boussole, color="blue") + 
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.25,8), ylim=c(43.2,43.75))


# For fish larvae
fl <- read.csv("visufront-fish-larvae.csv", sep=";", h=T)
head(fl)

sum(fl$larvae_nb) # 675 fish larvae
sum(station$eggs) # 548 eggs

# prepare the data
length(unique(fl$order)) # 13 orders
length(unique(fl$family)) # 29 families
species <- unique(str_c(str_c(fl$genera, fl$sp, sep=" "), str_c(" (", fl$descriptor, ")")))
length(species)  # 47 species
species <- sort(species)
#write.table(species, "species-only.csv", row.names=F, col.names=F)
species <- read.csv("species.csv", h=T, sep=";")
species$species <- as.character(species$species)

# add a column with the species name in the fl
fl$species <- str_c(str_c(fl$genera, fl$sp, sep=" "), str_c(" (", fl$descriptor, ")"))

# add a column to the stations to join it to the fish abundance
station$station <- str_c("station_", station$station_nb)
d <- join(fl, station, by="station")
d <- join(d, species)

abund <- ddply(d, .(station,habitat), summarize, sum=sum(larvae_nb), abundance=sum(larvae_nb)/(vol_out[1]-vol_in[1]) * 10000, order=length(unique(order)), family=length(unique(family)), species=length(unique(sp)), lat=lat[1], lon=lon[1], date=date[1])
head(abund)

# for fish larvae RAW ABUNDANCE


# pdf("raw_abund.pdf", width=10, height=7)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=sum, colour=date), data=abund) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, colour="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), colour="gray60", data=Boussole, vjust=2, hjust=1, size=4) +
	scale_size("Nb of larvae") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))


ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=sum), colour="gray20", data=abund) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, colour="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), colour="gray60", data=Boussole, vjust=2, hjust=1, size=4) +
	scale_size("Nb of larvae") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))


# For larval fish RELATIVE ABUNDANCE
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=abundance), colour="gray20", data=abund) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Relative Abund") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))
	


# For larval fish diversity NUMBER OF SPECIES
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=species), colour="gray20", data=abund) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Nb of species") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))


# For larval fish diversity NUMBER OF FAMILIES
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=family), colour="gray20", data=abund) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Nb of families") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))
	
	
# For larval fish diversity NUMBER OF ORDERS
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=order), colour="gray20", data=abund) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Nb of orders") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))
	
	
	
	
# For larval fish diversity OTHER THAT COASTAL (pelagic, mesopelagic and benthopelagic)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=abundance), colour="gray20", data=abund[-which(abund$habitat %in% "coastal"), ]) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Abundance offshore") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))
	

# For larval fish diversity OTHER THAT COASTAL (pelagic, mesopelagic and benthopelagic)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=abundance), colour="gray20", data=abund[which(abund$habitat %in% "mesopelagic"), ]) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Abundance Meso") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))
	
	
	

# For larval fish diversity OTHER THAT COASTAL (pelagic, mesopelagic and benthopelagic)
ggplot() + 
	geom_contour(aes(z=-z, x=x, y=y), colour="gray70", data=bathyDF) +
	geom_polygon(aes(x=lon, y=lat), data=coast, fill="grey60") +  
	geom_point(aes(x=lon, y=lat, size=order), colour="gray20", data=abund[which(abund$habitat %in% "coastal"), ]) + 
	geom_point(aes(x=lon, y=lat), data=Boussole, color="gray60") + 
	geom_text(aes(x=lon, y=lat, label="Boussole"), data=Boussole, color="gray60", vjust=2, hjust=1, size=4) +
	scale_size("Abundance coastal") +
	scale_x_continuous("Latitude") +
	scale_y_continuous("Longitude") +
	coord_map(xlim=c(7.1,8.1), ylim=c(43.2,43.75))
	
dev.off()		
	
	



#       PROJECT VISUFRONT
*Created by Robin Faillettaz on date 14/02/2014
UPMC - Laboratoire d'Océanographie de Villefranche-sur-Mer (LOV)*

_______________
Project to analyse the results of the ISIIS data from the VISUFRONT cruise, focusing on the cross-current transects. There are several parts:

- Process of the physical data
- Process of the biological data
- Analyse of combined biological and physical data

_______________


###LIBRARIES 

*lib_zooprocess.R*

	- Read pids
	- Read learning set
	- Extract prediction variables
	- Add probablity of being sorted in a specific category
	- Remove the probabilities of each image of being sorted in that specific category

*lib_process.R*

	- Path to data directory
	- Read ts data from the Tethys
	- Read CTD data from ISIIS
	- Detect casts
	- Compute distances (from start, Vlfr, or shore)

*lib_confusion.R*

	- Prepare confusion matrix
	- Plot confusion matrix as a heatmap
	- Compute statistics (recall, precision, etc.)
	- Generate plot to compare prediction quality

*lib_plot.R*

	- Add nice color gradients
	- Interpolation over time and distance
	- Smooth interpolation
	- Plot PCA results


###PROCESS PHYSICAL DATA

- Physical data are split into ADCP, CTD and drifters data
- ADCP, TS and GPS data are from the Tethys sensors
- Drifters were CODE drifters
- CTD data are from ISIIS sensor
- Only upcasts are relevant due to error around surface in downcasts
- Spatial PCA is done with CTD data interpolated over distance


####GPS AND TS DATA

*gps-process.R*

	- Extract GPS data from all available sources
	- Compare different GPS tracks
	- Compute mean bearing per 1 min, 15 sec and 1 sec
	- Save datasets of GPS + bearings at these time steps

*maps.R*

	- Generate geographical maps of the cruise

**ts-process.R*

	- Read TS files from the Tethys
	- Split into transects
	- Save them into separated files


####DRIFTERS

*drifters-plot.R*

	- Compute drifter speed
	- Generate a plot with drifters positions, ship track, and drifters interpolated position at time t


####ISIIS CTD DATA

*isiis-real_time.R*

	- Plot CTD data on realtime during transects

*isiis-process.R*

	- Read ISIIS CTD data
	- Clean up aberrant values
	- Cut by transect and save into separated files

*isiis-interp.R*

	- Read isiis CTD data per transect
	- Smooth data 
	- Save each interpolated transect into a new file

*plot-transects.R*

	- Plot transects and identify casts

*spatial-pca-isiis.R*

	- Read physical data
	- Compute anomalies from the mean profil of each varaible
	- Do a spatial PCA to better represent the location of the front
	- The front is where the variability is the highest

####ADCP

*adcp-checksum.R*

	- Check end of line of ADCP data

*adcp-explore.R*

	- Find number of correct/missing bearings

*adcp-process.R*

	- Fix bearings in ADCP data
	- Identify missing bearings
	- Compute bearing from the previous and next position
	- Assign the new bearing to missing ones



###PROCESS BIOLOGICAL DATA

- Process biological data from ISIIS and from plankton nets
- Count of each category are bined per 1 m and transformed per m-3
- Explore relation with environmental variables with regression trees
- Explore vertical migration

*plankton-nets-visufront.R*

	- Check diversity per sample
	- Check species distribution, split into coastal and pelagics

*dat1-analysis.R*

	- Results presented at OSM 2014
	- Process dat1 files (=outputs from Zooprocess)
	- Transform into abundances per cubic meters, bined per 1 m depth
	- Join with physical data for further analyses
	- Overlay biological and physical data
	- Integrate abundances per profiles
	- Compare larval fish abundances from PLANKTON NETS and ISIIS

*analysis-env_abund.R*

	- Check relation between abundances and environmental variables
	- Boosted Regression Trees between each category and environmental variable
	- Regression Trees ~
	- Explore vertical diel migrations (=vertical distribution between day and night)

# install/load necessary packages
library(devtools)
#install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
library(geoNEON)
library(geosphere)

# define local path to NEON l1 Woody Vegetation data 
woody_path <- "~/Documents/earth-analytics/data/SJER_2017/NEON_struct-woody-plant/NEON.D17.SJER.DP1.10098.001.2015-04.basic.20171016T160520Z/NEON.D17.SJER.DP1.10098.001.vst_mappingandtagging.basic.20171016T160520Z.csv"

# load in situ tree data 
woody <- read.csv(woody_path)

# remove all rows without stem distance & azimuth data
woody <- woody[complete.cases(woody$stemAzimuth) & 
               complete.cases(woody$stemDistance),]
                 

# use the geoNEON R package to get location information for each woody veg entry
woody_loc <- def.extr.geo.os(woody, 'namedLocation')


# to calculate precise geolocations in UTM for each mapped stem
# as described on page 8 of the NEON veg structure user guide:
# http://data.neonscience.org/api/v0/documents/NEON_vegStructure_userGuide_vA

# get easting/northing of reference point ID
ref_east <- as.numeric(woody_loc$api.easting)
ref_north <- as.numeric(woody_loc$api.northing)
theta <- woody_loc$stemAzimuth * pi

# calculate easting and northing for each plant
# add new columns to the woody veg data frame 
woody$easting <- ref_east + 
              woody_loc$stemDistance * sin(theta)
woody$northing <- ref_north + 
               woody_loc$stemDistance * cos(theta)



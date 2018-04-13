# load packages 
library(devtools)
library(geoNEON)
library(sp)           # CRS
library(swfscMisc)    # circle.polygon
library(rgdal)        # writeOGR
library(dplyr)


########################### SETUP ##################################
# set working directory to neon-veg in your local environment
setwd("~/github/neon-veg")

# define path to directory containing NEON l1 Woody Vegetation data 
main_path <- "~/Documents/earth-analytics/data/SJER_2017/NEON_struct-woody-plant/" #desktop
#main_path <- "/Users/victoriascholl/CU-Boulder/earthlab/NEON_data/2017_SJER/NEON_struct-woody-plant/" #laptop

# define path and filename of output shapefile to be written
shp_filename <- "output/test_sjer_polygons"

#####################################################################


# load local functions written in external R files. 
source("locate_woody_veg.R")
source("woody_df_to_shp.R")
source("merge_vst_tables.R")
source("get_vst_crs.R")
source("list_tiles_with_plants.R")
source("apply_area_threshold.R")
source("polygon_overlap.R")

# loop through folders of field data with different dates
dirs <- list.dirs(path = main_path )
dirs <-dirs[ grepl("NEON.D17.SJER", dirs) ]


first_loop <- 1 # loop counter
for (woody_path in dirs){
  
  # mapping and tagging table (contains stem locations)
  woody_mapping_path <- paste(woody_path, 
                              list.files(path = woody_path, 
                                         pattern = "mappingandtagging"), 
                              sep="/")
  
  # apparent individual table (contains height and crown diameter)
  woody_individual_path <- paste(woody_path, 
                                 list.files(path = woody_path, 
                                            pattern = "apparentindividual"), 
                                 sep="/")
  
  # load situ stem locations table; 
  # calculate mapped UTM locations of plants from distance/azimuth
  woody <- read.csv(woody_mapping_path)
  woody_utm <- locate_woody_veg(woody)
  
  # load "vst_apparentindividual" table
  woody_individual <- read.csv(woody_individual_path)
  
  # match mapped stems from "vst_mappingandtagging" with structure data 
  # from "vst_apparentindividual" based on individualID 
  woody_vst <- merge_vst_tables(woody_utm, woody_individual)
  
  # combine woody veg structure data to a single data frame 
  if (first_loop == 1){
    woody_all <- woody_vst
    woody_mapping_all <- woody_utm
    woody_individual_all <- woody_individual
    first_loop <- 0
    
  } else {
    woody_all <- rbind(woody_all, woody_vst)
    woody_mapping_all <- rbind(woody_mapping_all, woody_utm)
    woody_individual_all <- rbind(woody_individual_all, woody_individual)
  }
}


# list the 1km x 1km tiles containing field data
tiles <- list_tiles_with_plants(woody_all)


# create coordinate reference system object based on
# UTM zone info in the "vst_plotperyear" table
crs <- get_vst_crs(woody_path)


# create circular polygon for each stem based on maxCanopyDiameter
woody_polygons <- woody_df_to_shp(df=woody_all, 
                                  coord_ref=crs,
                                  shrink=1,
                                  num_sides = 8,
                                  shp_filename = shp_filename)


# Apply area threshold: 
# since the goal of this work is to create polygons that outline individual 
# plants to extract their spectra, we only want to keep polygons that are 
# large enough to extract "pure" spectra. The spatial resolution of the RGB 
# digital camera data is 25cm x 25cm. The area of the hyperspectral imagery 
# and gridded LiDAR-derived CHM Ecocystem Structure product have a coarser 
# spatial resolution, 1m^2^ or 100cm x 100cm. An area threshold was set here 
# to keep only polygons with an area greater than or equal to 4 HS / CHM pixels
# based on the maxCrownDiameter (in units of meters)

woody_thresh <- apply_area_threshold(woody_all,
                                     nPix=4)

# delete/merge/clip overlapping polygons
polygon_overlap(woody_thresh, nPix=4)

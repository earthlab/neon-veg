# load packages 
library(devtools)
library(geoNEON)
library(sp)           
library(swfscMisc)    
library(rgdal)        
library(dplyr)


########################### SETUP ##################################

# set working directory to neon-veg in local environment
setwd("~/github/neon-veg")

# path to directory containing NEON l1 Woody Vegetation data 
main_path <- "~/Documents/NEON/SJER/NEON_struct-woody-plant/" 

# path and filename of output shapefile to be written
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
source("get_poly.R")

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

# remove polygons with area < 4 hyperspectral pixels 
woody_thresh <- apply_area_threshold(woody_all,
                                     nPix=4)

# create circular polygon for each stem based on maxCanopyDiameter
woody_polygons <- woody_df_to_shp(df=woody_thresh, 
                                  coord_ref=crs,
                                  shrink=1,
                                  num_sides = 8,
                                  shp_filename = shp_filename)

# delete/merge/clip overlapping polygons
woody_final <- polygon_overlap(woody_polygons, nPix=4)

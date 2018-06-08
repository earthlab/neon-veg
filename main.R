# load packages 
library(devtools)
library(geoNEON)
library(sp)           
library(swfscMisc)    
library(rgdal)        
library(dplyr)
library(rgeos)
library(stringr)
library(ggplot2) 
library(tidyr)
library(dplyr)
library(purrr)
library(broom)


########################### SETUP ##################################

main_path <- "SJER/woody_veg" 

# specify output directory path and filename of output shapefile to be written
out_dir <- "output/"

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
source("check_create_dir.R")
source("make_species_table.R")
source("df_to_shp_points.R")
source("apply_height_threshold.R")
source("allometry_height_diam.R")

# create output directory if it does not exist 
check_create_dir(out_dir)

# loop through folders of field data with different dates
dirs <- list.dirs(path = main_path )
dirs <- dirs[grepl("NEON.D", dirs)]


first_loop <- 1 # loop counter
for (woody_path in dirs) {
  
  # mapping and tagging table (contains stem locations)
  woody_mapping_path <- paste(woody_path, 
                              list.files(path = woody_path, 
                                         pattern = "mappingandtagging"), 
                              sep = "/")
  
  # apparent individual table (contains height and crown diameter)
  woody_individual_path <- paste(woody_path, 
                                 list.files(path = woody_path, 
                                            pattern = "apparentindividual"), 
                                 sep = "/")
  
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
  if (first_loop == 1) {
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

# create text file to keep track of the number of trees after each step
count_file <- file(paste(out_dir,"tree_counts.txt", sep=""), "w")

# number of trees with mapped locations 
tree_count <- paste(as.character(nrow(woody_mapping_all)),
                    "trees with mapped locations",
                    sep=" ")
write(tree_count, count_file, append=TRUE)

# list the 1km x 1km tiles containing field data
tiles <- list_tiles_with_plants(woody_all, out_dir)

# create coordinate reference system object based on
# UTM zone info in the "vst_plotperyear" table
coord_ref <- get_vst_crs(woody_path)

# write shapefile with points for every mapped stem 
df_to_shp_points(woody_mapping_all, 
                 coord_ref, 
                 shp_filename = paste(out_dir,
                                      "mapped_stems",
                                      sep = ""))

# remove duplicate entries; keep most recent
woody_no_duplicates <- woody_all %>% 
  group_by(individualID) %>%
  slice(which.max(as.Date(date)))

# number of trees with complete entries, no duplicates 
# (location, species, height, crown diam)
tree_count <- paste(as.character(nrow(woody_no_duplicates)),
                    "trees with complete entries",
                    sep=" ")
write(tree_count, count_file, append=TRUE)

# write merged entries to csv
write.csv(woody_no_duplicates, file = paste(out_dir,"vst_merged.csv"))



# before applying area threshold, create polygon shapefile
# for all complete entries 
woody_df_to_shp(df = woody_no_duplicates, 
                     coord_ref = coord_ref,
                     shrink = 1,
                     num_sides = 8,
                     shp_filename = paste(out_dir,
                                           "polygons_all",
                                           sep = ""))

# remove polygons with area < 4 hyperspectral pixels 
woody_thresh <- apply_area_threshold(woody_no_duplicates,
                                     nPix = 4)

# number of trees after applying area threshold
tree_count <- paste(as.character(nrow(woody_thresh)),
                    "trees after applying area threshold",
                    sep=" ")
write(tree_count, count_file, append=TRUE)

# display table of species count in thresholded data frame 
species_table <- make_species_table(woody_thresh)
species_table

# create circular polygon for each stem based on max crown diameter
woody_polygons <- woody_df_to_shp(df = woody_thresh, 
                                  coord_ref = coord_ref,
                                  shrink = 1,
                                  num_sides = 24,
                                  shp_filename = paste(out_dir,
                                                       "polygons_filtered",
                                                       sep = ""))

# delete/merge/clip overlapping polygons
woody_final <- polygon_overlap(woody_polygons,
                               nPix = 4, 
                               shp_filename = paste(out_dir,
                                                    "polygons_checked_overlap",
                                                    sep = ""))

# write shapefile of mapped stem locations for final polygons 
stems_final <- as.data.frame(woody_final) 
idx_ID <- woody_thresh$individualID %in% woody_final$individualID
stems_final$easting <- woody_thresh$easting[idx_ID]
stems_final$northing <-  woody_thresh$northing[idx_ID]
df_to_shp_points(stems_final, 
                 coord_ref, 
                 shp_filename = paste(out_dir,
                                      "mapped_stems_final",
                                      sep = ""))


# number of trees after checking for overlap & 
# applying area threshold to clipped polygons 
tree_count <- paste(as.character(nrow(woody_final)),
                    "trees after checking for polygon overlap",
                    sep=" ")
write(tree_count, count_file, append=TRUE)

# close file keeping track of tree counts
close(count_file)

# construct height-crown diameter allometry for each species
allometries <- allometry_height_diam(stems_final)


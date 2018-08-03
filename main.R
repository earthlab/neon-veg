# load packages 
library(devtools)
library(geoNEON)
library(neonUtilities)
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

# specify the directory containing NEON data.
data_dir <- "NIWO/" 

# specify output directory path and filename of output shapefile to be written
out_dir <- "output/"

#####################################################################

# load local external functions 
source("supporting_functions.R")

# create output directory if it does not exist 
check_create_dir(out_dir)

# create text file to keep track of the number of trees after each step
count_file <- file(paste(out_dir,"tree_counts.txt", sep=""), "w")



# read & pre-process woody veg structure data -----------------------------


# define the path to the zipped woody veg data
woody_veg_filename = paste0(data_dir, 'woody_veg/NEON_struct-woody-plant.zip')

# use the stackByTable function to unzip the woody veg structure data and 
# combine the data into a single series of data tables. 
# "dpID" is the Data Product ID; use 'DP1.10098.001' for woody veg structure. 
neonUtilities::stackByTable(dpID = 'DP1.10098.001',
                             filepath = woody_veg_filename)


# read the mappingandtagging data table
woody_mapping_all = read.csv(paste0(data_dir,
                                    'woody_veg/',
                                    'NEON_struct-woody-plant/',
                                    'stackedFiles/',
                                    'vst_mappingandtagging.csv'))

# remove any duplicate individualID entries; keep most recent 
woody_mapping <- remove_duplicates(woody_mapping_all)

# check how many trees there are in the input mapping_and_tagging file 
tree_count <- paste0(as.character(length(unique(woody_mapping$individualID))), 
                     ' unique tree IDs in the mapping_and_tagging file')
print(tree_count)
write(tree_count, count_file, append=TRUE)


# calculate mapped UTM locations of plants from distance/azimuth
woody_utm <- locate_woody_veg(woody_mapping_all)

# check how many trees there are with mapped stem distance & azimuth
tree_count <- paste0(as.character(length(unique(woody_utm$individualID))), 
             ' trees with stem distance & azimuth values in the',
             ' mappingandtagging data table')
print(tree_count)
write(tree_count, count_file, append=TRUE)


# read the apparentindividual data table (contains height and crown diameter)
woody_individual_all = read.csv(paste0(data_dir,
                                    'woody_veg/',
                                    'NEON_struct-woody-plant/',
                                    'stackedFiles/',
                                    'vst_apparentindividual.csv'))

# remove any duplicate individalID entries; keep most recent 
woody_individual <- remove_duplicates(woody_individual_all)

# keep only the entries with crown diameter and tree height
# since these measurements are needed to create and compare polygons
woody_ind_complete <- woody_individual[complete.cases(woody_individual$height) & 
                                       complete.cases(woody_individual$maxCrownDiameter),]

# check how many trees there are in the input mapping_and_tagging file 
tree_count <- paste0(as.character(length(unique(woody_ind_complete$individualID))), 
             ' trees with height and max crown diameter values in the',
             ' apparentindividual data table')
print(tree_count)
write(tree_count, count_file, append=TRUE)

# match structural measurements from apparent_individual with 
# mapping_and_tagging entries based on individualID.
# remove any entries that have NA values for height or crown diameter. 
woody_merged <- merge(woody_utm, 
                      woody_ind_complete, 
                      by = "individualID") %>% 
  dplyr::filter(!is.na(height) & !is.na(maxCrownDiameter))

# check how many trees have mapped locations (utm coordinates), height, & diam
tree_count <- paste0(as.character(length(unique(woody_merged$individualID))), 
             ' trees with location, height, and max crown diameter values')
print(tree_count)
write(tree_count, count_file, append=TRUE)


# write merged mappingandtagging and apparentindividual entries to csv
write.csv(woody_merged, file = paste(out_dir,"vst_merged.csv"))

# list the 1km x 1km tiles containing mapped tree stems
tiles <- list_tiles_with_plants(woody_utm, out_dir)

# create coordinate reference system object based on
# UTM zone info in the "vst_plotperyear" data table
coord_ref <- get_vst_crs(paste0(data_dir,
                                'woody_veg/',
                                'NEON_struct-woody-plant/',
                                'stackedFiles/'))



# apply processing steps and create polygons ------------------------------

# write shapefile with points for every mapped stem 
df_to_shp_points(woody_utm, 
                 coord_ref, 
                 shp_filename = paste(out_dir,
                                      "mapped_stems",
                                      sep = ""))

# before applying area threshold, create polygon shapefile
# for all complete entries 
woody_df_to_shp(df = woody_merged, 
                     coord_ref = coord_ref,
                     shrink = 1,
                     num_sides = 24,
                     shp_filename = paste(out_dir,
                                           "polygons_all",
                                           sep = ""))

# remove polygons with area < 4 hyperspectral pixels 
woody_thresh <- apply_area_threshold(woody_merged,
                                     nPix = 4)

# number of trees after applying area threshold
tree_count <- paste0(as.character(length(unique(woody_thresh$individualID))),
                    " trees after applying area threshold")
print(tree_count)
write(tree_count, count_file, append=TRUE)

# create circular polygon for each stem based on max crown diameter
woody_polygons <- woody_df_to_shp(df = woody_thresh, 
                                  coord_ref = coord_ref,
                                  shrink = 1,
                                  num_sides = 24,
                                  shp_filename = paste(out_dir,
                                                       "polygons_filtered",
                                                       sep = ""))

# delete/clip overlapping polygons
woody_final <- polygon_overlap(woody_polygons,
                               nPix = 4, 
                               shp_filename = paste(out_dir,
                                                    "polygons_checked_overlap",
                                                    sep = ""))

# write shapefile of mapped stem locations for final polygons 
# get the easting and northing tree coordinates from woody_thresh data frame
stems_final <- as.data.frame(woody_final) 
stems_coordinates <- woody_thresh %>% dplyr::select(individualID, easting, northing)
stems_final <- merge(stems_final, stems_coordinates, by = "individualID")
df_to_shp_points(stems_final, 
                 coord_ref, 
                 shp_filename = paste(out_dir,
                                      "mapped_stems_final",
                                      sep = ""))


# number of trees after checking for overlap & 
# applying area threshold to clipped polygons 
tree_count <- paste0(as.character(length(unique(woody_final$individualID))),
                    " trees after checking for polygon overlap")
print(tree_count)
write(tree_count, count_file, append=TRUE)

# close file keeping track of tree counts
close(count_file)

# display table of species count in the merged complete data
species_table <- make_species_table(stems_final)
species_table

# construct height-crown diameter allometry for each species
allometries <- allometry_height_diam(stems_final)

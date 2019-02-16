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

# path to NEON Woody Plant Vegetation Structure zip file 
data_path <- "data/NIWO/NEON_struct-woody-plant.zip" 

# specify output directory path and filename of output shapefile to be written
out_dir <- "output/"

# factor to specify the relative size of the polygons generated.
# the radius of each polygon is divided by this value. 
# when this factor == 1, the polygons have a radius of (maxCrownDiameter / 2),
# so the polygons are created with radii half the size as the maxCrownDiameter. 
# when this factor == 2, the polygons have a radius of (maxCrownDiameter / 4),
# so the polygons are one quarter of the size as the maxCrownDiameter.  
crown_size_factor <- 1

# number of pixels used to threshold the area of polygons.
# polygons smaller than this will be excluded at certain steps during the analysis. 
# the units of this threshold are # of pixels in the hyperspectral data set (1m x 1m) 
area_thresh_pixels <- 4

#####################################################################

# load local external functions 
source("supporting_functions.R")

# create output directory if it does not exist 
check_create_dir(out_dir)

# create text file to keep track of the number of trees after each step
count_file <- file(paste(out_dir,"tree_counts.txt", sep=""), "w")



# read & pre-process woody veg structure data -----------------------------


# use the stackByTable function to unzip the woody veg structure data and 
# combine the data into a single series of data tables. 
# "dpID" is the Data Product ID; use 'DP1.10098.001' for woody veg structure. 
neonUtilities::stackByTable(dpID = 'DP1.10098.001',
                            filepath = data_path)


# read the mappingandtagging data table
woody_mapping_all = read.csv(paste0(tools::file_path_sans_ext(data_path),
                                    '/stackedFiles/',
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
woody_individual_all = read.csv(paste0(tools::file_path_sans_ext(data_path),
                                       '/stackedFiles/',
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
coord_ref <- get_vst_crs(paste0(tools::file_path_sans_ext(data_path),
                                '/stackedFiles/'))


# remove multi-bole entries -----------------------------------------------

# identify multi-bole trees: sometimes, there are multiple entries with identical
# coordinates, height, and crown diameter (but their individualID's are different,
# with "A", "B", etc. appended on the end of the last five numbers.) In this step, 
# the individualID values are assessed to find multi-bole sets. 
# if the coordinates, height, and crown diameters are identical, then the
# entry with a letter appended on the end is deleted from the analysis. This prevents
# duplicate polygons being used to extract spectra. 

# get a list of all individual ID strings 
individualIDs_all <- as.character(unique(droplevels(woody_merged$individualID)))

# split each individual ID using a period delimiter
last_id_section <- sapply(stringr::str_split(individualIDs_all, "[.]"), tail, 1)

# create another list without any letters
last_digits <- gsub("[^0-9]","",last_id_section)

# create a boolean vector where bole entries (which contain letters)
# are True, but stem entries (without letters) are False 
is_bole <- last_id_section != last_digits

# create a lookup table with the id information
id_lut <- as.data.frame(last_digits) %>% 
            mutate(individualIDs_all = individualIDs_all,
                   last_id_section = last_id_section,
                   is_bole = is_bole, 
                   height = woody_merged$height,
                   maxCrownDiameter = woody_merged$maxCrownDiameter)

# count the frequency of each ID number. identify the ID's with more than one entry. 
multiple_ids <- as.data.frame(table(last_digits)) %>% 
  filter(Freq >1)

# create a list to populate with individualIDs to remove
remove_ids <- c()

# loop through the ID's that appear more than once in the data set
for(id in as.character(multiple_ids$last_digits)){
  
  print(id)
  
  # get the complete individual IDs 
  duplicates <- print(id_lut[id_lut$last_digits == id,])
  
  # see if the height and diameter values are identical 
  if(var(duplicates$height)==0 && var(duplicates$maxCrownDiameter) == 0){
    remove_ids <- c(remove_ids, duplicates$individualIDs_all[duplicates$is_bole==TRUE])
  }
  
}

# remove the entries with the multi-bole individualIDs identified in the previous step. 
woody_multibole_removed <- woody_merged %>% 
  filter(!(individualID %in%remove_ids))

# check how many trees are left after removing multi-bole entries
tree_count <- paste0(as.character(length(unique(woody_multibole_removed$individualID))), 
                     ' trees remaining after multi-bole entries were removed')
print(tree_count)
write(tree_count, count_file, append=TRUE)


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
                     shrink = crown_size_factor,
                     num_sides = 24,
                     shp_filename = paste(out_dir,
                                           "polygons_all",
                                           sep = ""))

# write polygons to file after multi-bole entries have been removed 
# (this is before the area threshold filtering step)
woody_df_to_shp(df = woody_multibole_removed, 
                coord_ref = coord_ref,
                shrink = crown_size_factor,
                num_sides = 24,
                shp_filename = paste(out_dir,
                                     "polygons_multibole_removed",
                                     sep = ""))


# remove polygons with area < 4 hyperspectral pixels 
woody_thresh <- apply_area_threshold(woody_multibole_removed,
                                     nPix = area_thresh_pixels)

# number of trees after applying area threshold
tree_count <- paste0(as.character(length(unique(woody_thresh$individualID))),
                    " trees after applying area threshold")
print(tree_count)
write(tree_count, count_file, append=TRUE)

# create circular polygon for each stem based on max crown diameter
# after polygons have been filtered based on area threshold
woody_polygons <- woody_df_to_shp(df = woody_thresh, 
                                  coord_ref = coord_ref,
                                  shrink = crown_size_factor,
                                  num_sides = 24,
                                  shp_filename = paste(out_dir,
                                                       "polygons_filtered",
                                                       sep = ""))

# remove engulfed polygons 
woody_delete_engulfed <- delete_engulfed(woody_polygons,
                                         paste(out_dir, 
                                               "polygons_deleted_engulfed"))

# count how many polygons are left after removing the engulfed ones. 
tree_count <- paste0(as.character(dim(woody_delete_engulfed)[1]),
                     " trees after deleting engulfed polygons")
print(tree_count)
write(tree_count, count_file, append=TRUE)



# Clip/delete overlapping polygons ----------------------------------------

# clip overlapping polygons (delete them if remaining polygons < area threshold) 
woody_clipped <- clip_overlap(woody_delete_engulfed,
                              nPix = area_thresh_pixels,
                              shp_filename = paste0(out_dir,
                                                   "polygons_clipped_overlap"))


# write shapefile of mapped stem locations for final polygons 
# get the easting and northing tree coordinates from woody_thresh data frame
stems_final <- as.data.frame(woody_clipped) 
stems_coordinates <- woody_thresh %>% dplyr::select(individualID, easting, northing)
stems_final <- merge(stems_final, stems_coordinates, by = "individualID")
df_to_shp_points(stems_final, 
                 coord_ref, 
                 shp_filename = paste(out_dir,
                                      "mapped_stems_final",
                                      sep = ""))


# number of trees after checking for overlap & 
# applying area threshold to clipped polygons 
tree_count <- paste0(as.character(length(unique(woody_clipped$individualID))),
                    " trees after checking for polygon overlap")
print(tree_count)
write(tree_count, count_file, append=TRUE)

# close file keeping track of tree counts
close(count_file)

# display table of species count in the merged complete data
species_table <- make_species_table(stems_final)
species_table

# construct height-crown diameter allometry for each species.
allometries <- allometry_height_diam(stems_final, 
                                     paste0(out_dir,"allometry.png"))



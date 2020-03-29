# This script extracts features (remote sensing data) for each
# sample (pixel) within the specified shapefile (containing points or
# polygons that correspond to trees at the NEON site)




# loop through stacked AOP data and extract features 
stacked_aop_list <- list.files(file.path(dir_data_out, "stacked_aop_data"),
                               full.names = TRUE)

# clip data cube - extract features ----------------------------------------

print("Extracting features for tree points / polygons in:  ") 
print(shapefile_filename)

# read the shapefile layer 
shp <- sf::st_read(shapefile_filename)

# add columns for the center location of each tree 
shp_coords <- shp %>% 
  sf::st_centroid() %>%  # get the centroids first for polygon geometries 
  sf::st_coordinates() %>% 
  as.data.frame()

# add new columns for the tree location coordinates 
shp$X <- shp_coords$X
shp$Y <- shp_coords$Y

# get a description of the shapefile to use for naming outputs
shapefile_description <- tools::file_path_sans_ext(basename(shapefile_filename))

# loop through AOP tiles 
for (stacked_aop_filename in stacked_aop_list) {
  
  print("Extracting features in tile:")
  print(basename(stacked_aop_filename))
  
  # read current tile of stacked AOP data 
  stacked_aop_data <- readRDS(stacked_aop_filename)
  
  # construct the easting northing string for naming outputs
  east_north_string <- paste0(stacked_aop_data$eastingIDs[1],
                              "_",
                              stacked_aop_data$northingIDs[1])
  
  # figure out which trees are within the current tile by comparing each
  # X,Y coordinate to the extent of the current tile 
  trees_in <- shp %>% 
    dplyr::filter(X >= extent(stacked_aop_data)[1] & 
                    X < extent(stacked_aop_data)[2] & 
                    Y >= extent(stacked_aop_data)[3] & 
                    Y  < extent(stacked_aop_data)[4])
  
  print(paste0(as.character(nrow(trees_in))," trees in current tile"))
  
  # if no polygons are within the current tile, skip to the next one
  if (nrow(trees_in)==0){
    print("no trees located within current tile... skipping to next shapefile")
    next
  }
  
  # convert from SF obect to Spatial object for clipping
  trees_in_sp <- sf::as_Spatial(trees_in,
                                IDs = as.character(trees_in$indvdID))
  
  # clip the hyperspectral raster stack with the polygons within current tile.
  # the returned objects are data frames, each row corresponds to a pixel in the
  # hyperspectral imagery. The ID number refers to which tree that the 
  # the pixel belongs to. A large polygon will lead to many extracted pixels
  # (many rows in the output data frame), whereas tree stem points will
  # lead to a single extracted pixel per tree. 
  print("Extracting features for each tree from the data cube... ")
  extracted_spectra <- raster::extract(stacked_aop_data, 
                                       trees_in_sp, 
                                       df = TRUE)
  
  # VS-NOTE TO DO: 
  # adjust this extract step to only get pixels WITHIN each tree polygon,
  # also try calculating the percentage that each pixel is within a polygon
  # and keep only pixels with > 50% overlap 
  
  # merge the extracted spectra and other data values with the tree info 
  tree_metadata <- data.frame(trees_in) %>% 
    mutate(ID = 1:nrow(trees_in))
  
  # create a list of increasing integer counts to keep track of how many rows 
  # (pixels or spectra) belong to each tree 
  for (j in unique(extracted_spectra$ID)){
    if(j==1){
      counts = 1:sum(extracted_spectra$ID==j)
    }
    else{
      counts = append(counts, 1:sum(extracted_spectra$ID==j))
    }
  }
  
  # combine the additional data with each spectrum for writing to file.
  # remove the geometry column to avoid issues when writing to csv later 
  spectra_write <- merge(tree_metadata,
                         extracted_spectra,
                         by="ID") %>% 
    dplyr::mutate(spectra_count = counts)%>% 
    dplyr::select(ID, spectra_count, everything()) %>% 
    dplyr::select(-geometry)
  
  # write extracted spectra and other remote sensing data values to file 
  write.csv(spectra_write, 
            file = file.path(dir_data_out,
                          paste0("extracted_features_",
                          east_north_string, "_",
                          shapefile_description,
                          ".csv"))) 
}  

# combine all extracted features into a single .csv
csvs <- list.files(dir_data_out, full.names = TRUE)

# refine the output csv selection 
csvs <- csvs[grepl(paste0("*000_", shapefile_description, ".csv"), csvs)]

# combine all .csv data into a single data frame 
for (c in 1:length(csvs)){
  print(csvs[c])
  csv <- read.csv(csvs[c])
  
  if(c==1){
    spectra_all <- csv
  } else {
    # add a bias value to the ID column, so in the end
    # the ID values will range from 1 to n_trees
    csv$ID <- csv$ID + max(spectra_all$ID)
    spectra_all <- rbind(spectra_all, csv)
  }
}

# remove the unneccessary column "X.1"
spectra_all <- spectra_all %>% dplyr::select(-X.1)

# write ALL the spectra to a single .csv file 
# VS-NOTE return this filename for use in subsequent scripts? 
extracted_features_filename <- file.path(dir_data_out,
                                  paste0(shapefile_description,
                                         "-extracted_features.csv"))
write.csv(spectra_all,
          file=extracted_features_filename)

# delete the individual csv files for each tile 
file.remove(csvs)

  
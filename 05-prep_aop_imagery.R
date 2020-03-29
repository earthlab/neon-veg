# This script stacks up the AOP remote sensing imagery and prepares it for
# extracting descriptive features for classification. 


# Specify the paths for each data directory
h5_dir <- file.path(dir_data_raw, "hyperspectral")
chm_dir <- file.path(dir_data_raw, "chm")
slope_dir <- file.path(dir_data_raw, "slope")
aspect_dir <-file.path(dir_data_raw, "aspect")
rgb_dir <- file.path(dir_data_raw, "rgb")
veg_indices_dir <- file.path(dir_data_raw, "veg_indices")

# output directory for stacked AOP data
stacked_aop_data_dir <- file.path(dir_data_out, "stacked_aop_data")
check_create_dir(stacked_aop_data_dir)


# list the files in each data directory; filter results based on file type
# hyperspectral data - list the .h5 files 
h5_list <- list.files(path = h5_dir, full.names = TRUE) 
chm_list <- list.files(path = chm_dir, full.names = TRUE)
slope_list <- list.files(path = slope_dir, full.names = TRUE)
aspect_list <- list.files(path = aspect_dir, full.names = TRUE)
rgb_list <- list.files(path = rgb_dir, full.names = TRUE)
veg_indices_list <- list.files(path = veg_indices_dir, full.names = TRUE)
# remove any entries from the list that are .zip files
veg_indices_list<- veg_indices_list[grepl("*.tif$", veg_indices_list)]
veg_indices_names <- c("ARVI","EVI","NDLI","NDNI","NDVI","PRI","SAVI")



# create data cubes with AOP-derived features 

start_time <- Sys.time() # start the timer 

# loop through the tiles; build up a data cube for each tile 

for (h5 in h5_list) {
  
  # print current hyperspectral filename 
  print(h5)
  
  # each current hyperspectral tile must be read and stacked into a 
  # georeferenced rasterstack object (so it can be clipped with point / polygon
  # shapefiles). The process of creating a rasterstack takes a while for 
  # each tile, so after creating each rasterstack once, each object gets 
  # written to a file. 
  
  # Build up the rasterstack filename by parsing out the easting/northing
  # coordinates from the current h5 filename.
  # (this rasterstack just contains the hyperspectral layers)
  
  # parse the UTM easting and northing values from the current h5 filename
  easting <- stringr::str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][5]
  northing <- stringr::str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][6]
  # combine them with an underscore; use this to find corresponding tiles 
  # of various remote sensing data
  east_north_string <- paste0(easting,"_",northing)
  
  # generate a filename for the stacked AOP data
  stacked_aop_data_filename = file.path(stacked_aop_data_dir,
                              paste0("stacked_aop_data_",
                                     east_north_string, ".rds"))
  # check if a .rds file already exists for the current feature data cube
  if (file.exists(stacked_aop_data_filename)){
    
    # if it exists, read that instead of re-generating the same rasterstack.
    message("stacked_aop_data already created for current tile.")
    
    # restore / read the rasterstack from file
    #stacked_aop_data <- readRDS(file = stacked_aop_data_filename)
    
  } else{
    
    # if it doesn't exist, create the features from the aop data to file 
    message("creating stacked_aop_data for current tile...")
    
    
    
    # hyperspectral and lidar features ----------------------------------------
    
    # Build up the h5 rasterstack filename
    rasterstack_filename <- paste0(h5_dir, "rasterstack_",
                                   east_north_string, ".rds")
    
    print(paste("rasterstack filename: ", rasterstack_filename))
    
    # check to see if a .rds file already exists for the current tile
    # within the hyperspectral directory - just HS reflectance bands so far 
    if (file.exists(rasterstack_filename)){
      
      # if it exists, read that instead of re-generating the same rasterstack.
      message("reading hyperspectral rasterstack (already created for current tile)...")
      # restore / read the rasterstack from file
      s <- readRDS(file = rasterstack_filename)
      
    } else{
      
      # if it doesn't exist, generate the rasterstack. 
      message("creating hyperspectral rasterstack for current tile...")
      # create a georeferenced rasterstack using the current hyperspectral tile
      s <- stack_hyperspectral(h5, out_dir = dir_data_out)
      # save the rasterstack to file 
      #saveRDS(s, file = rasterstack_filename)
    }
    
    # read the corresponding remote sensing data layers for current tile
    print("Reading lidar-derived layers for the training data cube...")
    chm <- raster::raster(grep(east_north_string, chm_list, value=TRUE))
    slope <- raster::raster(grep(east_north_string, slope_list, value=TRUE))
    aspect <- raster::raster(grep(east_north_string, aspect_list, value=TRUE))
    
    # for the vegetation indices, go into the corresponding folder for current tile
    # and get a list of all the vegIndex geotiffs. then read all of those geotiffs 
    # into a single raster stack.
    print("Reading vegetation indices for the training data cube...")
    veg_indices <- raster::stack(grep(east_north_string, 
                                      veg_indices_list, 
                                      value=TRUE))
    
    # set the raster name for each layer to be simply the name of the data 
    # (i.e. "aspect") as opposed to the full filename 
    # (i.e. ""NEON_D13_NIWO_DP3_452000_4431000_aspect")
    names(chm) <- "chm"
    names(slope) <- "slope"
    names(aspect) <- "aspect"
    # name each of the vegetation index layers based on the last piece of each 
    # respective filename, e.g. "NDVI" and 
    names(veg_indices) <- sapply(stringr::str_split(names(veg_indices),"_"),tail,1)
    
    # RGB features ------------------------------------------------------------
    
    # The RGB data tile has 10,000 x 10,000 pixels, 10cm spatial resolution. 
    # All other layers tiles have 1,000 x 1,000 pixels, 1 meter spatial resolution. 
    # Aggregate red, green, blue intensity within each coarser grid cell using
    # statistics such as mean and standard deviation. 
    
    rgb_features_filename <- file.path(rgb_dir, paste0(
                              "rgb_features_",east_north_string, ".rds"))
    
    print(paste("rgb_features: ", rgb_features_filename))
    
    # check if a .rds file already exists for the current feature data cube
    if (file.exists(rgb_features_filename)){
      
      # if it exists, read that instead of re-generating the same rasterstack.
      message("reading rgb_features (already created for current tile)...")
      
      # restore / read the rasterstack from file
      rgb_features <- readRDS(file = rgb_features_filename)
      
    } else{
      
      # if it doesn't exist, create the features from the aop data to file 
      message("creating rgb_features for current tile...")
      
      # rgb data has 3 bands. read each one individually 
      rgb_red <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), band = 1) 
      rgb_green <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), band = 2) 
      rgb_blue <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), band = 3) 
      
      # The "fact" parameter of the raster::aggregate function is the number of cells
      # in each direction (horizontal and vertically) to aggregate across.
      # Since the RGB data has a spatial resolution that is 1/10th of the 
      # other data layers (10cm compared to 1m), fact should be 10 to produce 
      # an output raster with 1000 x 1000 pixels. 
      # mean intensity per 1m x 1m grid cell 
      rgb_meanR <- raster::aggregate(rgb_red, fact = 10, fun = mean)
      rgb_meanG <- raster::aggregate(rgb_green, fact = 10, fun = mean)
      rgb_meanB <- raster::aggregate(rgb_blue, fact = 10, fun = mean)
      
      # standard deviation of intensity per 1m x 1m grid cell 
      rgb_sdR <- raster::aggregate(rgb_red, fact = 10, fun = sd)
      rgb_sdG <- raster::aggregate(rgb_green, fact = 10, fun = sd)
      rgb_sdB <- raster::aggregate(rgb_blue, fact = 10, fun = sd)
      
      # (mean + SD) gives some idea about relative position on the spectrum
      # and variation rather than just variation (standard deviation on its own)
      rgb_mean_sd_R <- (rgb_meanR + rgb_sdR)
      rgb_mean_sd_G <- (rgb_meanG + rgb_sdG)
      rgb_mean_sd_B <- (rgb_meanB + rgb_sdB)
      
      
      # future work: would be awesome to calculate texture features: 
      # https://cran.r-project.org/web/packages/radiomics/vignettes/TextureAnalysis.html
      
      
      # to confirm the order of the red, green, and blue intensities,
      # stack all 3 bands into a RasterStack and use the plotRGB function.
      # the colors appear natural, so r=1, g=2, b=3
      #rgb_stack <- raster::stack(rgb_red,rgb_green, rgb_blue)
      #raster::plotRGB(rgb_stack, r = 1, g = 2, b = 3)
      
      # set the names of each layer to reflect the metric it contains
      names(rgb_meanR) <- "rgb_meanR" # mean intensity 
      names(rgb_meanG) <- "rgb_meanG"
      names(rgb_meanB) <- "rgb_meanB"
      names(rgb_sdR) <- "rgb_sdR"     # standard deviation of intensity 
      names(rgb_sdG) <- "rgb_sdG"
      names(rgb_sdB) <- "rgb_sdB"
      names(rgb_mean_sd_R) <- "rgb_mean_sd_R" # mean plus standard deviation
      names(rgb_mean_sd_G) <- "rgb_mean_sd_G"
      names(rgb_mean_sd_B) <- "rgb_mean_sd_B"
      
      # stack up all the RGB features
      rgb_features <- raster::stack(rgb_meanR, rgb_meanG, rgb_meanB,
                                    rgb_sdR, rgb_sdG, rgb_sdB,
                                    rgb_mean_sd_R, rgb_mean_sd_G, rgb_mean_sd_B)
      
      # save the stacked features from the aop data to file 
      #saveRDS(rgb_features, file = rgb_features_filename)
    }
    
    # Create the pixel number grid as a layer to add to the data cube. 
    # this one keeps track of individual pixel ID's
    # to avoid duplicate spectra being extracted. Basically, assign an integer ID
    # to each pixel in the 1000x1000 raster. This raster needs to have the same 
    # dimensions, extent, crs as the other layers so they can be stacked together. 
    # create a vector of IDs from 1 to the number of pixels in one band (#rows x #cols)
    pixelID <- 1:(nrow(s) * ncol(s))
    # add tile east, north coordinates - how to do this if raster values must be numeric? 
    #pixelID <- paste(pixelID, east_north_string, sep="_") 
    # reshape this 1D vector into a 2D matrix 
    dim(pixelID) <- c(nrow(s),ncol(s))
    # create a raster layer of pixel numbers 
    pixelNumbers <- raster::raster(pixelID, crs = crs(s))
    extent(pixelNumbers) <- extent(s)
    names(pixelNumbers) <- "pixelNumber"
    
    # Create similar layers to keep track of the tile where the pixel is located
    print("Creating layers for easting and northing per pixel...")
    eastingID <- rep(as.numeric(easting), times = (nrow(s) * ncol(s)))
    northingID <- rep(as.numeric(northing), times = (nrow(s) * ncol(s)))
    # reshape to be two-dimensional
    dim(eastingID) <- c(nrow(s),ncol(s))
    dim(northingID) <- c(nrow(s),ncol(s))
    # create rasters to contain the easting and northing values
    eastingIDs <- raster::raster(eastingID, crs = crs(s))
    northingIDs <- raster::raster(northingID, crs = crs(s))
    # assign extent and CRS to match the other layers in the stack
    extent(eastingIDs) <- extent(s)
    extent(northingIDs) <- extent(s)
    names(eastingIDs) <- "eastingIDs"
    names(northingIDs) <- "northingIDs"
    
    
    # now, all of the remote sensing data files have been read in for the current
    # tile. add each one to the hyperspectral data stack along with the 
    # layer to keep track pixel number within the tile. 
    stacked_aop_data <- raster::addLayer(s, chm, slope, aspect, veg_indices, 
                                         rgb_features, pixelNumbers, 
                                         eastingIDs, northingIDs)
    print("Stacked AOP data for current tile. ")
    
    # save the stacked AOP data to file for easy clipping later
    saveRDS(stacked_aop_data, file = stacked_aop_data_filename)
    
  }

}

end_time <- Sys.time()
elapsed <- end_time - start_time
print("Elapsed time: ")
print(elapsed)

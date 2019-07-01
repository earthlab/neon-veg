# classify_species.R 
# Repo: https://github.com/earthlab/neon-veg 
# 
# This script quantifies species classification accuracy using 
# freely-available NEON data and Random Forest classification. 
# 
# First, download Woody Vegetation Structure data for your NEON site and 
# and time period of interest. Run the "create_tree_points_polygons.R" script
# to create points and polygons for each individual tree in the data set. 
# 
# Before running this script, download the following NEON airborne remote 
# sensing mosaic products, from which predictive features will be extracted 
# for the species classification. Place each type of imagery in its own
# folder within the "data" folder, named accordingly: 
#   hyperspectral (Spectrometer orthorectified surface directional reflectance)
#   chm (Ecosystem structure, Canopy Height Model (CHM))
#   rgb (High-resolution orthorectified camera imagery mosaic)
#   aspect (Slope and Aspect - LiDAR)
#   slope (Slope and Aspect - LiDAR)
#   vegIndices (Vegetation indices - spectrometer - mosaic)
#
#
# Author: Victoria Scholl
# Last updated: 06/30/19


# Installation & setup ----------------------------------------------------

# load necessary R packages 
library(rhdf5)
library(rgdal)
library(raster)
library(ggplot2)
library(tidyr)
library(sf)
library(dplyr)
library(data.table) 
library(stringr) # str_split function
library(randomForest) # randomForest function
library(ggbiplot) # for PCA visualization with a biplot
library(rfUtilities) # to calculate OA, PA, UA accuracies 

# set working directory
setwd("~/github/neon_veg/")

# load local functions in external files 
source("code/supporting_functions.R")

# code for NEON site 
site_code <- "NIWO"

# define the "bad bands" wavelength ranges in nanometers, where atmospheric 
# absorption creates unreliable reflectance values. 
bad_band_window_1 <- c(1340, 1445)
bad_band_window_2 <- c(1790, 1955)

# taxon ID's to predict
taxonList <- c("ABLAL","PICOL","PIEN","PIFL2")

# define the output directory. If it doesn't exist already, create it.
check_create_dir("output/") # create top level "output" directory
out_dir <- paste0("output/", site_code, "/")
check_create_dir(out_dir) # create output folder for site
check_create_dir(paste0(out_dir, "figures/")) # create output figures folder


# shapefiles to test: tree points, polygons -------------------------------

# Shapefile set(s) to associate a tree species labels for each pixel. 
# These shapefiles were generated using the "create_tree_points_polygons.R"
# workflow within this neon_veg repo, where each entry (point or polygon) 
# corresponds to a tree in the woody vegetation in-situ field NEON data set.

# directory with shapefiles to test (tree stem locations and crown polygons)
shapefile_dir <- paste0("output/", site_code, "/")

# define the subdirectory (destination or dsn when reading shapefile)
# and the layer name (the filename before the .shp extension) 
# in a vector for each of the shapefile scenarios: 

# (1) "allStems" - all stem points at NIWO that also have crown diameter 
# measurements. this includes multi-bole entries, which result in 
# duplicated spectra. 
allStems_layer <- c("allStems",
                    "shapefiles_maxDiameter/",
                    "mapped_stems_with_crown_diameter")

# (2) "all_polygonspolygons with max crown diameter, one generated for each 
# stem point at NIWO. Just as for the "allStems" layer above, the shapfile 
# used here includes multi-bole entries, since this is how the data come when
# downloaded straight from the data portal. 
allPolygons_maxDiameter_layer <- c("allPolygons_maxDiameter",
                                   "shapefiles_maxDiameter/",
                                   "polygons_all")

# (3) "allPolygons_halfDiameter" - polygons with half the maximum diameter, 
# one generated for each stem point at NIWO. 
allPolygons_halfDiameter_layer <- c("allPolygons_halfDiameter",
                                    "shapefiles_halfDiameter/",
                                    "polygons_all")

# (4) "neonvegStems_maxDiameter" - stem point locations after being processed 
# by the neon_veg workflow steps. These stem locations correspond to polygons 
# generated with maximum diameter. 
neonvegStems_maxDiameter_layer <- c("neonvegStems_maxDiameter",
                                    "shapefiles_maxDiameter/",
                                    "mapped_stems_final")

# (5) "neonvegPolygons_halfDiameter" -  neon_veg workflow polygons generated 
# with half the maximum crown diameter per individual tree. 
neonvegPolygons_halfDiameter_layer <- c("neonvegPolygons_halfDiameter",
                                        "shapefiles_halfDiameter/",
                                        "polygons_clipped_overlap")

# (6) "neonvegPolygons_maxDiameter" - neon_veg workflow polygons generated 
# with maximum crown diameter per individual tree. 
neonvegPolygons_maxDiameter_layer <- c("neonvegPolygons_maxDiameter",
                                       "shapefiles_maxDiameter/",
                                       "polygons_clipped_overlap")

# create a data frame containing the description, directory, and shapefile name
# for each of the shapefile scenarios to be tested 
shapefile_layer_names <- as.data.frame(rbind(allStems_layer,
                                           allPolygons_halfDiameter_layer,
                                           allPolygons_maxDiameter_layer,
                                           neonvegStems_maxDiameter_layer,
                                           neonvegPolygons_halfDiameter_layer,
                                           neonvegPolygons_maxDiameter_layer),
                                     stringsAsFactors = FALSE) %>% 
  `colnames<-`(c("description", "dsn", "layer")) 
rownames(shapefile_layer_names) <- 1:nrow(shapefile_layer_names)

print("shapefile layers to test:")
shapefile_layer_names



# VS-NOTE: add a vector where the user can specify which features to extract



# Input AOP airborne remote sensing data --------------------------------------

# Assemble all of the NEON AOP remote sensing data layers into a giant data 
# cube. A set of spectral and structural features from this data cube will 
# be extracted for each sample (pixel) that intersects with the shapefile 
# entries. These data were downloaded for the site of interest from 
# https://data.neonscience.org
# VS-NOTE: In the future, access the NEON API to automatically download these.

# Specify the paths for each data directory:
# hyperspectral .h5 and .rds
h5_dir <- paste0("data/", site_code, "/hyperspectral/") 
# CHM geotiffs 
chm_dir <- paste0("data/", site_code, "/chm/")      
# slope geotiffs
slope_dir <- paste0("data/", site_code, "/slope/")   
# aspect geotiffs
aspect_dir <- paste0("data/", site_code, "/aspect/") 
# rgb image geotiffs
rgb_dir <- paste0("data/", site_code, "/rgb/")
# vegetation index .tifs 
vegIndices_dir <- paste0("data/", site_code, "/vegIndices/") 
# output directory for stacked AOP data
stacked_aop_data_dir <- paste0("data/",site_code, "/stacked_aop_data/")
check_create_dir(stacked_aop_data_dir)

# list the files in each data directory; filter results based on file type
# hyperspectral data - list the .h5 files 
h5_list <- list.files(path = h5_dir, full.names = TRUE)
h5_list <- h5_list[grepl("*.h5", h5_list)] 
# canopy height model data - list the .tif files 
chm_list <- list.files(path = chm_dir, full.names = TRUE)
chm_list <- chm_list[grepl("*CHM.tif$", chm_list)]
# lidar-derived slope data - list the .tif files 
slope_list <- list.files(path = slope_dir, full.names = TRUE)
slope_list <- slope_list[grepl("*slope.tif$", slope_list)]
# lidar-derived aspect data - list the .tif files
aspect_list <- list.files(path = aspect_dir, full.names = TRUE)
aspect_list <- aspect_list[grepl("*aspect.tif$", aspect_list)]
# digital camera rgb data - list the .tif files 
rgb_list <- list.files(path = rgb_dir, full.names = TRUE)
rgb_list<- rgb_list[grepl("*image.tif$", rgb_list)]
# hyperspectral-derived veg indices - list the subdirectories for each index 
vegIndices_list <- list.dirs(path = vegIndices_dir, full.names = TRUE)
vegIndices_names <- c("ARVI","EVI","NDLI","NDNI","NDVI","PRI","SAVI")

# read wavelengths from text file. If this file does not exist,
# create this text file using the first hyperspectral image file. 
wavelengths_filename <- paste0(out_dir,"wavelengths.txt")
if (file.exists(wavelengths_filename)) {
  wavelengths = as.numeric(unlist(read.table(wavelengths_filename,
                                             sep="\n",
                                             skip = 1,
                                             col.names = 'wavelength')))
} else {
  # get the name of the first hyperspectral file 
  h5 <- h5_list[1]
  # list the contents of HDF5 file
  h5_struct <- rhdf5::h5ls(h5, all=T)
  # read the wavelengths of the hyperspectral image bands
  wavelength_tag <- paste0(h5_struct$group[grepl(
    "/Reflectance/Metadata/Spectral_Data", h5_struct$group)][1],"/Wavelength")
  wavelengths <- rhdf5::h5read(h5,wavelength_tag)
  # write the exact wavelengths to file for future use 
  write.table(data.frame(wavelengths = wavelengths),
              wavelengths_filename,
              sep="\n",
              row.names=FALSE)
}




# Read and stack AOP remote sensing data layers ------------------------------

start_time <- Sys.time() # start the timer 

# loop through the tiles; build up a data cube for each tile 
for (h5 in h5_list) {
  
  # print current hyperspectral filename 
  print(h5)
  
  # Each current hyperspectral tile must be read and stacked into a 
  # georeferenced rasterstack object (so it can be clipped with point / polygon
  # shapefiles). The process of creating a rasterstack takes a while for 
  # each tile, so after creating each rasterstack once, each object gets 
  # written to a file. 
  
  # Build up the rasterstack filename by parsing out the easting/northing
  # coordinates from the current h5 filename.
  # (this rasterstack just contains the hyperspectral layers)
  
  # parse the UTM easting and northing values from the current h5 filename
  # based on this format: "NEON_D13_NIWO_DP3_451000_4432000_reflectance.h5"
  easting <- stringr::str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][5]
  northing <- stringr::str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][6]
  
  # combine them with an underscore; use this to find corresponding tiles 
  # of various remote sensing data
  east_north_string <- paste0(easting, "_", northing)
  
  stacked_aop_data_filename = paste0(stacked_aop_data_dir, "stacked_aop_data_",
                                     east_north_string, ".rds")
  # check if a .rds file already exists for the current feature data cube
  if (file.exists(stacked_aop_data_filename)){
    
    # if it exists, read that instead of re-generating the same rasterstack.
    message("reading stacked_aop_data (already created for current tile)...")
    
    # restore / read the rasterstack from file
    stacked_aop_data <- readRDS(file = stacked_aop_data_filename)
    
  } else{
    
    
    # Hyperspectral layers --------------------------------------------------
    
    # Build up the h5 rasterstack filename 
    rasterstack_filename <- paste0(h5_dir, "rasterstack_",
                                   east_north_string, ".rds")
    
    # check to see if a .rds file already exists for the current tile
    # within the hyperspectral directory - just HS reflectance bands so far 
    if (file.exists(rasterstack_filename)){
      
      # if it exists, read that instead of re-generating the same rasterstack.
      message("reading rasterstack (already created for current tile)...")
      # restore / read the rasterstack from file
      s <- readRDS(file = rasterstack_filename)
      
    } else{
      
      # generate the rasterstack
      message("creating rasterstack for current tile...")
      
      # create a georeferenced rasterstack using the current hyperspectral tile
      s <- stack_hyperspectral(h5, out_dir)
      
      # save the rasterstack to file 
      saveRDS(s, file = rasterstack_filename)
      
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
    vegIndices <- raster::stack(list.files(grep(east_north_string, 
                                                vegIndices_list, 
                                                value=TRUE), 
                                           pattern="tif$", full.names=TRUE))
    
    # set the raster name for each layer to be simply the name of the data 
    # (i.e. "aspect") as opposed to the full filename 
    # (i.e. ""NEON_D13_NIWO_DP3_452000_4431000_aspect")
    names(chm) <- "chm"
    names(slope) <- "slope"
    names(aspect) <- "aspect"
    # name each of the vegetation index layers based on the last piece of each 
    # respective filename, e.g. "NDVI" and 
    names(vegIndices) <- sapply(stringr::str_split(names(vegIndices), 
                                                   "_"), tail, 1)
    
    # High-res RGB data -----------------------------------------------------
    
    # The RGB data tile has 10,000 x 10,000 pixels, 10cm spatial resolution. 
    # All other layers tiles have 1,000 x 1,000 pixels, 1 meter spatial resolution. 
    # Aggregate red, green, blue intensity within each coarser grid cell using
    # statistics such as mean and standard deviation. 
    
    rgb_features_filename <- paste0(rgb_dir, "rgb_features_",
                                    east_north_string, ".rds")
    
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
      rgb_red <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), 
                                band = 1) 
      rgb_green <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), 
                                  band = 2) 
      rgb_blue <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), 
                                 band = 3) 
      
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
      
      
      # VS-NOTE: in the future, it would be awesome to calculate texture features: 
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
      saveRDS(rgb_features, file = rgb_features_filename)
    }
    
    # Create the pixel number grid as a layer to add to the data cube. 
    # this one keeps track of individual pixel ID's
    # to avoid duplicate spectra being extracted. Basically, assign an integer ID
    # to each pixel in the 1000x1000 raster. This raster needs to have the same 
    # dimensions, extent, crs as the other layers so they can be stacked together. 
    # create a vector of IDs from 1 to the number of pixels in one band (#rows x #cols)
    pixelID <- 1:(nrow(s) * ncol(s))
    # Add tile east, north coordinates - how to do this if raster values must be numeric? 
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
    
    
    # Now, all of the remote sensing data files have been read in for the current
    # tile. add each one to the hyperspectral data stack along with the 
    # layer to keep track pixel number within the tile. 
    stacked_aop_data <- raster::addLayer(s, chm, slope, aspect, vegIndices, 
                                         rgb_features, pixelNumbers, 
                                         eastingIDs, northingIDs)
    print("Stacked remote sensing layers for current tile. ")
    
    # save the stacked AOP data to file for easy clipping later
    saveRDS(stacked_aop_data, file = stacked_aop_data_filename)
    
  }
  
  
  
  
  
  
  
  
  # plot rasters ------------------------------------------------------------
  # show raster examples for specific tile with a road (defining feature)
  if(east_north_string == "452000_4432000"){
    
    # plot the high-res RGB image after it has been aggregated to coarser grid
    rgb_stack <- raster::stack(stacked_aop_data$rgb_meanR, 
                               stacked_aop_data$rgb_meanG, 
                               stacked_aop_data$rgb_meanB)
    rgb_brick <- raster::brick(rgb_stack)
    raster::plotRGB(rgb_brick,
                    r = 1, g = 2, b = 3,
                    stretch = "lin",
                    axes = TRUE,
                    main = paste("Digital Camera RGB Composite \nfor tile",
                                 easting, northing),
                    cex.main = 1)
    
    # plot RGB composite using selected bands from the hyperspectral imagery
    r <- s[[which.min(abs(wavelengths - 620))]] # 620nm Red
    g <- s[[which.min(abs(wavelengths - 555))]] # 555nm Green
    b <- s[[which.min(abs(wavelengths - 450))]] # 450nm Blue 
    hs_rgb_stack <- raster::stack(r,g,b)
    hs_rgb_brick <- raster::brick(hs_rgb_stack)
    
    raster::plotRGB(hs_rgb_brick,
                    r = 1, g = 2, b = 3,
                    stretch = "lin",
                    axes = TRUE,
                    main = paste("Hyperspectral-derived \nRGB Composite using",
                                 "bands 620nm, 555nm, 450nm\n for tile",
                                 easting, northing),
                    #xlab = "Easting (m)",
                    #ylab = "Northing (m)", 
                    cex.main = 1)
    
    # write the RGB composite to an image file for visualizing in QGIS
    writeRaster(hs_rgb_brick,
                paste0(out_dir, "figures/", "rgb_composite_", 
                       east_north_string, ".tif"), 
                format = "GTiff",
                overwrite = TRUE)
    
    # CHM 
    raster::plot(chm,
                 col = grey(1:100/100),
                 axes = TRUE,
                 main = paste("Lidar-Derived Canopy Height Model \nfor tile",
                        easting, northing),
                 xlab = "Easting (m)",
                 ylab = "Northing (m)",
                 cex.main = 1, # title text size 
                 # format the legend on the color bar
                 legend.args = list(text = "Height above ground [m]",side = 2, 
                                    font=2, line = 0.5, cex=0.8))
    
    # slope
    raster::plot(slope,
                 col = grey(1:100/100),
                 axes = TRUE,
                 main = paste("Lidar-Derived Slope\nfor tile",
                              easting, northing),
                 xlab = "Easting (m)",
                 ylab = "Northing (m)",
                 cex.main = 1,
                 legend.args = list(text="Ratio of rise over run [degrees]",
                                  side=2, font=2, 
                                  line=0.5, cex=0.8)) # legend on color bar)
    
    # aspect 
    raster::plot(aspect,
                 col=grey(1:100/100),
                 axes = TRUE,
                 main = paste("Lidar-Derived Aspect \nfor tile",
                              easting, northing),
                 xlab = "Easting (m)",
                 ylab = "Northing (m)",
                 cex.main = 1,
                 legend.args = list(text = paste("Direction of steepest", 
                                                 "slope \n [degrees from",
                                                 "North]"),
                                               side = 2, font = 2, 
                                  line = 0.5, cex = 0.8)) # legend on color bar)
    
    # NDVI vegetation index
    raster::plot(vegIndices[["NDVI"]],
                 col = grey(1:100/100),
                 axes = TRUE,
                 main = paste("Hyperspectral-Derived NDVI \nfor tile",
                              easting, northing),
                 xlab = "Easting (m)",
                 ylab = "Northing (m)",
                 cex.main = 1,
                 legend.args = list(text = "NDVI [0,1]",side = 2, font = 2, 
                                    line = 0.5, cex = 0.8))
  }
  
  
  # Clip data cube - extract features -------------------------------------
  
  # Loop through each of the shapefile sets 
  for(i in 1:nrow(shapefile_layer_names)){ 
    print(paste0("Currently extracting features for tree points / ",
                 "polygons in:  ", 
                 shapefile_layer_names$description[i]))
    
    # read the shapefile layer 
    shp <- rgdal::readOGR(dsn = paste0(shapefile_dir,shapefile_layer_names$dsn[i]),
                          layer = shapefile_layer_names$layer[i])
    
    # convert to SF object
    shp_sf <- sf::st_as_sf(shp)
    
    # add columns for the center location of each tree 
    shp_coords <- shp_sf %>% 
      sf::st_centroid() %>%  # get the centroids first for polygon geometries 
      sf::st_coordinates() %>% 
      as.data.frame()
    
    # add new columns for the tree location coordinates 
    shp_sf$X <- shp_coords$X
    shp_sf$Y <- shp_coords$Y
    
    # figure out which trees are within the current tile by comparing each
    # X,Y coordinate to the extent of the current tile 
    trees_in <- shp_sf %>% 
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
    
    # VS-NOTE: TO DO: 
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
    spectra_write <- base::merge(tree_metadata,
                           extracted_spectra,
                           by="ID") %>% 
      dplyr::mutate(spectra_count = counts)%>% 
      dplyr::select(ID, spectra_count, everything()) %>% 
      dplyr::select(-geometry)
    
    # write extracted spectra and other remote sensing data values to file 
    write.csv(spectra_write, 
              file = paste0(out_dir,
                            "extracted_features_",
                            east_north_string, "_",
                            shapefile_layer_names$description[i], ".csv")) 
  }
  
}

end_time <- Sys.time()
elapsed <- end_time - start_time
print("Elapsed time: ")
print(elapsed)




#############
############








# combine all extracted features into a single 
# .csv for each of the shapefile sets.
for(i in 1:nrow(shapefile_layer_names)){
  
  csvs <- list.files(out_dir, full.names = TRUE)
  
  # refine the output csv selection 
  csvs <- csvs[grepl(paste0("*000_", shapefile_layer_names$description[i], ".csv"), csvs)]
  
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
  write.csv(spectra_all,
            file=paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                        shapefile_layer_names$description[i],".csv"))
  
}















# ribbon plots ------------------------------------------------------------

# loop through each shapefile name, read the .csv containing spectral reflectance 
# for all trees within each data set, generate a ribbon plot and write to image file
# in the figures/ output directory 
for(i in 1:nrow(shapefile_layer_names)){
  extracted_features_filename <- paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                                        shapefile_layer_names$description[i],".csv")
  
  createRibbonPlot(wavelengths, extracted_features_filename)
  
  createSeparateRibbonPlots(wavelengths, extracted_features_filename)
  
}


















# count number of polygons and pixels per shapefile scenario ------------------
for(i in 1:nrow(shapefile_layer_names)){
  print(i)
  extracted_features_filename <- paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                                        shapefile_layer_names$description[i],".csv")
  print("Currently counting the number of individual IDs and extracted spectra in:")
  print(extracted_features_filename)
  
  # read the values extracted from the data cube
  df_orig <- read.csv(extracted_features_filename)
  
  # filter down to the species of interest 
  df_orig <- df_orig %>% 
    dplyr::filter(taxonID %in% taxonList)
  
  # for the first iteration, create a row of values 
  if(i ==1){
    count <- c(length(unique(df_orig$indvdID)), length(unique(df_orig$pixelNumber)))
    
    # append to the initial row of count values in subsequent iterations 
  } else{
    count <- rbind(count, c(length(unique(df_orig$indvdID)), length(unique(df_orig$pixelNumber))))
  }
  
  
  # for the comparison between classification accuracies where each
  # pair of shapefiles (the "raw NEON" version compared to the corresopnding
  # "neon_veg" filtered version), create variables to keep track of the 
  # individual tree ID's that exist within each neon_veg shapefile. 
  # check for the neon_veg stems shapefile; keep track of the individual ID's
  if(grepl("neonvegStems", extracted_features_filename)){
    print(paste0("Recording individualIDs for neonveg Stems:", extracted_features_filename))
    neonveg_stem_IDs <- unique(droplevels(df_orig$indvdID))
  }
  # check for the neon_veg halfDiameter shapefile; keep track of the individual ID's 
  if(grepl("neonvegPolygons", extracted_features_filename) & 
     grepl("halfDiam", extracted_features_filename)){
    print(paste0("Recording individualIDs for neonveg halfDiam:", extracted_features_filename))
    neonveg_halfDiam_IDs <- unique(droplevels(df_orig$indvdID))
  }
  # check for the neon_veg maxDiameter shapefile; keep track of the individual ID's 
  if(grepl("neonvegPolygons", extracted_features_filename) & 
     grepl("maxDiam", extracted_features_filename)){
    print(paste0("Recording individualIDs for neonveg maxDiam:", extracted_features_filename))
    neonveg_maxDiam_IDs <- unique(droplevels(df_orig$indvdID))
  }
}
countDF <- data.frame(count,
                      row.names = shapefile_layer_names$description)
colnames(countDF) <- c("nIndvdID", "nPixelNumbers")









###########################################################################
# Random Forest Classification --------------------------------------------
###########################################################################


# At this point, all of the shapefile scenarios have been used to extract
# features from the giant remote sensing data cube.
outDescription <- "rf_allSamplesPerClass_ntree5000_allBandRefl_nVar6_mean-sd-RGB_independentValidationSet20percent/"
outDescription <- "rf_neonvegIDsForBothShapefiles_ntree5000_pca2InsteadOfWavelengths_nVar6_mean-sd-RGB_independentValidationSet20percent/" 
outDescription <- "rf_allSamplesPerClass_ntree5000_pca2InsteadOfWavelengths_nVar6_mean-sd-RGB_independentValidationSet20percent/" 

check_create_dir(paste0(out_dir,outDescription))

# RF tuning parameter, number of trees to grow. deafualt value 500
ntree <- 5000

# boolean variable. if TRUE, select random minSamples per class to reduce bias
randomMinSamples <- FALSE

# to remove the sample size bias, if TRUE this filters down each of the raw NEON 
# shapefile data sets to only contain the individualIDs present in the neon_veg set 
neonvegIDsForBothShapefiles <- FALSE

# boolean variable. if TRUE, keep separate set for validation
independentValidationSet <- TRUE 
# randomly select this amount of data for training, use the rest for validation
percentTrain <- 0.8 

pcaInsteadOfWavelengths <- FALSE
nPCs <- 2 # number of PCAs to keep 

# keep most important variables and run RF again with reduced feature set 
#keepMostImpVar <- FALSE

# create boxplots and write to file 
createBoxplots <- TRUE

# open a text file to record the output results 
rf_output_file <- file(paste0(out_dir,outDescription,
                              "rf_model_summaries.txt"), "w")

# create an empty matrix to summarise the model Accuracies
rfAccuracies <- data.frame(matrix(ncol = 4, nrow = nrow(shapefile_layer_names)))
colnames(rfAccuracies) <- c("shapefileDescription", "OA_OOB", "OA_IndVal","K")
rfAccuracies$shapefileDescription <- shapefile_layer_names$description

# create an empty matrix to summarise the variable importance rankings 
nVar = 6 # the n most important variables to record for each shapefile
rfVarImp <- data.frame(matrix(ncol = 1 + (2*nVar), nrow = nrow(shapefile_layer_names)))
colnames(rfVarImp) <- c("shapefileDescription", paste0("MDA_", as.character(1:nVar)),
                        paste0("MDG_", as.character(1:nVar)))
rfVarImp$shapefileDescription <- shapefile_layer_names$description


# start the timer
start_time <- Sys.time()

# read wavelengths if not previously created
wavelengths = as.numeric(unlist(read.table(paste0(out_dir,"wavelengths.txt"),
                                           sep="\n",
                                           skip = 1,
                                           col.names = "wavelength")))


if(independentValidationSet){
  # randomly select <percentTrain> of data 
  # from the neon_veg half diameter polygons for training,
  # use the remaining samples for validation.
  # keep track of which pixelNumber, easting, northing
  # that these pixels correspond to. 
  
  # read the file with spectra to sample randomly for validation trees 
  validationSourceFilename <- paste0(out_dir, site_code, 
                                     "_spectral_reflectance_ALL_",
                                     shapefile_layer_names$description[5],
                                     ".csv")
  # remove any spectra with a height of 0
  # and remove any factors
  df_val <- read.csv(validationSourceFilename) %>% 
    dplyr::filter(chm>0) %>% 
    base::droplevels() 
  
  # combine the pixelNumber, easting, and northing into a string
  # for each row with underscore separation
  dfIDs <- paste(df_val[,"pixelNumber"],
                 df_val[,"eastingIDs"],
                 df_val[,"northingIDs"], 
                 sep = "_") 
  
  # add this as a new column to the data frame 
  df_val <- df_val %>% mutate(dfIDs = dfIDs)
  
  # remove any duplicated spectra from consideration
  df_val <- df_val[!duplicated(df_val$dfIDs),]
  
  # randomly sample rows from this data set 
  set.seed(104)
  # the number of sampled rows is calculated based on 
  # percentTrain and the number of rows in the validation set. 
  # percentTrain may have a value like 0.80 (80% data used to train)
  train <- sample(nrow(df_val), 
                  percentTrain*nrow(df_val), 
                  replace = FALSE)
  # identify the other 0.20 (20%) of samples for independent validation
  validationSet <- df_val[-train,]
  
  # combine the pixelNumber, easting, and northing into a string
  # for each row with underscore separation
  valIDs <- validationSet$dfIDs
  
  # in the loop below, when each set of training features is 
  # prepared for the classifier, any spectra with
  # pixelNumer_eastingIDs_northingIDs that belong in the validation set
  # will be removed. 
  
} else{
  # otherwise, just remove the pixelNumber, eastingIDs and northingIDs
  # since they should not be input to the classifier. 
  features <- features %>% 
    dplyr::select(-c(pixelNumber, eastingIDs, northingIDs))
}


# loop through all shapefile sets 
for(i in 1:nrow(shapefile_layer_names)){ 
  
  # filename of current .csv file to read 
  extracted_features_filename <- paste0(out_dir, site_code, 
                                        "_spectral_reflectance_ALL_",
                                        shapefile_layer_names$description[i],
                                        ".csv")
  print("Currently training the random forest with features extracted using:")
  print(extracted_features_filename)
  
  # read the values extracted from the data cube
  df_orig <- read.csv(extracted_features_filename)
  
  # filter down to the species of interest 
  df_orig <- df_orig %>% 
    dplyr::filter(taxonID %in% taxonList)
  
  
  # testing the influence of sampling bias 
  if(neonvegIDsForBothShapefiles){
    # filter the raw NEON data down to contain the same individual ID's as the 
    # corresponding neon_veg data set. This will remove the potential influence
    # of sampling bias (since there is nearly double # of samples for the raw
    # NEON data sets) and just compare any influence that the clipping/filtering
    # steps have on the polygons 
    if(grepl("allStems", extracted_features_filename)){
      print(paste0("Filtering individualIDs to match those in neonveg_Stems:", extracted_features_filename))
      filterIDs <- neonveg_stem_IDs
    }
    # check for the neon_veg halfDiameter shapefile; keep track of the individual ID's 
    if(grepl("allPolygons", extracted_features_filename) & 
       grepl("halfDiam", extracted_features_filename)){
      print(paste0("Filtering individualIDs to match those in neonveg_halfDiam:", extracted_features_filename))
      filterIDs <- neonveg_halfDiam_IDs
    }
    # check for the neon_veg maxDiameter shapefile; keep track of the individual ID's 
    if(grepl("allPolygons", extracted_features_filename) & 
       grepl("maxDiam", extracted_features_filename)){
      print(paste0("Filtering individualIDs to match those in neonveg maxDiam:", extracted_features_filename))
      filterIDs <- neonveg_maxDiam_IDs
    }
    df_orig <- df_orig %>% dplyr::filter(indvdID %in% filterIDs) %>% droplevels()
  }
  
  
  
  # Remove any spectra that have a height == 0
  print(paste0(as.character(sum(df_orig$chm==0)), 
               " pixels have a height of 0 in the CHM"))
  print("Removing these rows from the training set ... ")
  
  # also reset the factor levels (in case there are dropped taxonID levels)
  df <- df_orig %>% filter(chm>0) %>% droplevels()
  
  # remove the bad bands from the list of wavelengths 
  remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                                 wavelengths < bad_band_window_1[2]) | 
                                (wavelengths > bad_band_window_2[1] & 
                                   wavelengths < bad_band_window_2[2])]
  
  # create a LUT that matches actual wavelength values with the column names,
  # X followed by the rounded wavelength values. Remove the rows that are 
  # within thebad band ranges. 
  wavelength_lut <- data.frame(wavelength = wavelengths,
                               xwavelength = paste0("X", 
                                                    as.character(round(wavelengths))),
                               stringsAsFactors = FALSE) %>% 
    filter(!wavelength %in% remove_bands) 
  
  # features to use in the RF models.
  # this list is used to filter the columns of the data frame,
  # to remove the ones containing other metadata per tree from the model. 
  featureNames <- c("taxonID", 
                    wavelength_lut$xwavelength,
                    "chm", 
                    "slope", 
                    "aspect",
                    "ARVI",
                    "EVI",
                    "NDLI",
                    "NDNI",
                    "NDVI",
                    "PRI",
                    "SAVI",
                    #"rgb_meanR",
                    #"rgb_meanG",
                    #"rgb_meanB",
                    "rgb_mean_sd_R",
                    "rgb_mean_sd_G",
                    "rgb_mean_sd_B",
                    "pixelNumber",
                    "eastingIDs",
                    "northingIDs")
  
  # filter the data to contain only the features of interest 
  features <- df %>% 
    dplyr::select(featureNames)
  
  # testing whether PCA yields better accuracy than individual wavelength reflectance data
  if(pcaInsteadOfWavelengths == TRUE){
    
    # remove the individual spectral reflectance bands from the training data
    features <- features %>% dplyr::select(-c(wavelength_lut$xwavelength))
    
    print(colnames(features))
    
    # PCA: calculate Principal Components 
    hs <- df %>% dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()
    hs_pca <- stats::prcomp(hs, center = TRUE, scale. = TRUE)
    summary(hs_pca)
    features <- cbind(features, hs_pca$x[,1:nPCs]) # add first n PCs to features data frame
    # visualize where each sample falls on a plot with PC2 vs PC1 
    ggbiplot::ggbiplot(hs_pca, 
                       choices = 1:2, # which PCs to plot
                       obs.scale = 1, var.scale = 1, # scale observations & variables 
                       var.axes=FALSE, # remove arrows 
                       groups = df$taxonID, # color the points by species
                       ellipse = TRUE, # draw ellipse around each group
                       circle = TRUE # draw circle around center of data set
    )   + 
      ggtitle("PCA biplot, PC1 and PC2") + 
      scale_color_brewer(palette="Spectral") + 
      theme_bw()
    # save to file 
    ggplot2::ggsave(paste0(out_dir, outDescription, 
                           "pcaPlot_",shapefile_layer_names$description[i],
                           ".png"))
    
  }
  
  print("Features used in current RF model: ")
  print(colnames(features))
  
  
  # count the number of samples per species 
  featureSummary <- features %>%
    dplyr::group_by(as.character(taxonID)) %>%
    dplyr::summarize(total = n()) 
  
  print("number of samples per species class")
  print(featureSummary)
  
  # randomly select <percentTrain> of data for training,
  # use the remaining samples for validation
  if(independentValidationSet){
    
    if(i == 1){
      
      # filter the data to contain only the features of interest 
      validationSet <- validationSet %>% 
        dplyr::select(featureNames)
      
      if(pcaInsteadOfWavelengths){
        # perform PCA
        # remove the individual band reflectances
        wl_removed <- validationSet %>% as.data.frame() %>% 
          dplyr::select(-c(wavelength_lut$xwavelength))
        
        # isolate the individual band reflectances
        val_hs <- validationSet %>% 
          dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()
        
        val_hs_pca <- stats::prcomp(val_hs, center = TRUE, scale. = TRUE)
        summary(val_hs_pca)
        validationSet <- cbind(wl_removed, val_hs_pca$x[,1:nPCs])
      }
      
    }
    
    # concatenate pixel #, easting, and northing 
    # for an individual entry per pixel in the input imagery,
    # i.e. "264231_452000_4431000"
    features$dfIDs <- paste(features[,"pixelNumber"],
                            features[,"eastingIDs"],
                            features[,"northingIDs"], 
                            sep = "_") 
    
    print("The following pixelNumber_easting_northing values are ")
    print("found in the current input set and the validation set: ")
    print(intersect(features$dfIDs, valIDs))
    
    print("Originally this many rows in the features DF: ")
    print(nrow(features))
    
    # remove any spectra that are being kept separate for the 
    # independent validation set 
    features <- features[!(features$dfIDs %in% valIDs),]
    
    print("After removing pixels that are in the validation set, nrow = ")
    print(nrow(features))
    
    # remove the pixelNumber, easting, and northing columns since they
    # are not input features to the train the classifier 
    features <- features %>% 
      dplyr::select(-c(pixelNumber, eastingIDs, northingIDs, dfIDs))
    
  }
  
  if(randomMinSamples){
    # reduce number of samples per species to avoid classifier bias
    
    # count the minimum number of samples for a single class
    minSamples <- min(featureSummary$total)  
    print(paste0("Randomly selecting ",
                 as.character(minSamples),
                 " samples per species class to avoid classifier bias"))
    
    # isolate the samples per species
    taxon1 <- features[features$taxonID==taxonList[1],]
    taxon2 <- features[features$taxonID==taxonList[2],]
    taxon3 <- features[features$taxonID==taxonList[3],]
    taxon4 <- features[features$taxonID==taxonList[4],]
    
    # keep random minSamples of each species; merge
    taxon1 <- taxon1[sample(nrow(taxon1), minSamples), ]
    taxon2 <- taxon2[sample(nrow(taxon2), minSamples), ]
    taxon3 <- taxon3[sample(nrow(taxon3), minSamples), ]
    taxon4 <- taxon4[sample(nrow(taxon4), minSamples), ]
    
    features <- rbind(taxon1, taxon2, taxon3, taxon4)
    
  } else{
    print("Using all samples per class")
  }
  
  
  
  # TRAIN RF CLASSIFIER using training set  ---------------------------------
  rf_startTime <- Sys.time()
  set.seed(104)
  rf_model <- randomForest::randomForest(as.factor(features$taxonID) ~ .,
                                         data=features, 
                                         importance=TRUE, 
                                         # ntree = # of trees to grow
                                         ntree=ntree) 
  print("randomForest time elapsed for model training: ")
  print(Sys.time()-rf_startTime)
  
  print(rf_model)
  
  # save RF model to file 
  save(rf_model, file = paste0(out_dir, outDescription,"rf_model_",
                               shapefile_layer_names$description[i],".RData"))
  
  # write all relevant information to the textfile: 
  # shapefile name
  write(shapefile_layer_names$description[i], rf_output_file, append=TRUE)
  write("\n", rf_output_file, append=TRUE) #newline
  
  # number of samples per class
  featureSummary <- data.frame(featureSummary)
  colnames(featureSummary) <- c("taxonID","numberOfSamples")
  capture.output(featureSummary, file = rf_output_file, append=TRUE)
  
  # y = predicted data; (horizontal axis)
  # x = observed data (true class labels) (vertical axis)
  accuracy <- rfUtilities::accuracy(x = rf_model$y,
                                    y = rf_model$predicted)
  
  # record each accuracy metric in the table for a final comparison.
  # round each value to the nearest decimal place 
  rfAccuracies$OA_OOB[i] <- round(accuracy$PCC, 1) # Overall Accuracy
  rfAccuracies$K[i] <- round(accuracy$kappa, 3) #Cohen's Kappa 
  
  # record the users and producer's accuracies for each specie s
  rfUsersProducers <- data.frame(matrix(ncol = 3, nrow = nrow(shapefile_layer_names)))
  #colnames(rfAccuracies) <- c("shapefileDescription", "OA", "K")
  #rfAccuracies$shapefileDescription <- shapefile_layer_names$description
  
  write("\nOverall Accuracy:", rf_output_file, append=TRUE) #newline
  write(accuracy$PCC, rf_output_file, append=TRUE)
  
  write("\nUser's Accuracy:", rf_output_file, append=TRUE) #newline
  capture.output(accuracy$users.accuracy, file = rf_output_file, append=TRUE)
  
  write("\nProducer's Accuracy:", rf_output_file, append=TRUE) #newline
  capture.output(accuracy$producers.accuracy, file = rf_output_file, append=TRUE)
  
  write("\nConfusion Matrix:", rf_output_file, append=TRUE) #newline
  capture.output(accuracy$confusion, file = rf_output_file, append=TRUE)
  
  write("\nCohen's Kappa:", rf_output_file, append=TRUE) #newline
  capture.output(accuracy$kappa, file = rf_output_file, append=TRUE)
  
  
  # INDEPENDENT VALIDATION  -------------------------------------------------
  
  # predict species ID for validation set 
  if(independentValidationSet){
    predValidation <- predict(rf_model, validationSet, type = "class")
    confusionTable <- table(predValidation, validationSet$taxonID)
    print(confusionTable)
    val_OA <- sum(predValidation == validationSet$taxonID) / 
      length(validationSet$taxonID)
    print(paste0("overall accuracy predicting validation set: ",
                 as.character(val_OA)))
    # write the accuracy summary data frame to file 
    write.csv(confusionTable,
              paste0(out_dir, outDescription, 
                     "rfConfusionMatrix_independentValidationSet_",
                     shapefile_layer_names$description[i],"_Accuracy_",
                     as.character(round(val_OA, 3)),".csv"))
    # write the Overall Accuracy for the independent validation set to the
    # Accuracies data frame
    rfAccuracies$OA_IndVal[i] <- round(val_OA*100, 1) # Overall Accuracy
  }
  
  
  # write all relevant information to the textfile: 
  # features used to describe each sample (pixel)
  write("\ndescriptive features used to train this model: ", rf_output_file, append=TRUE) #newline
  write(colnames(features), rf_output_file, append=TRUE)
  
  # RF model summary, OOB error rate 
  capture.output(rf_model, file = rf_output_file, append=TRUE)
  
  
  # VARIABLE IMPORTANCE  ----------------------------------------------------
  
  # What variables were important? --> Consult the variable importance plot. 
  varImportance <- data.frame(randomForest::importance(rf_model))
  varImportance$feature <- rownames(varImportance)
  varImportance <- varImportance %>% 
    select(feature, MeanDecreaseAccuracy, MeanDecreaseGini, everything())
  varImportanceMDA <- varImportance %>% dplyr::arrange(desc(MeanDecreaseAccuracy))
  varImportanceMDG <- varImportance %>% dplyr::arrange(desc(MeanDecreaseGini))
  
  # create bar plot to illustrate variable importance MDA
  ggplot(data = varImportanceMDA, aes(x = reorder(feature, MeanDecreaseAccuracy), 
                                      y = MeanDecreaseAccuracy,
                                      fill = MeanDecreaseAccuracy)) + 
    geom_bar(stat = 'identity', color = "black", size = 0.1, show.legend = FALSE) + 
    labs(x = "Variable", y = "Importance (MDA)") +
    coord_flip() + 
    theme_bw() + 
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9,"BuGn")),
                         values = rescale(varImportanceMDA$MeanDecreaseAccuracy, 
                                          to = c(0, 1))) + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) + 
    ggtitle("Variable Importance: Mean Decrease in Accuracy")
  
  # save variable importance bar plot to image file 
  ggsave(filename = paste0(out_dir, outDescription,"varImpPlot_MDA_",
                           shapefile_layer_names$description[i],
                           ".png"))
  
  
  # create bar plot to illustrate variable importance MDG
  ggplot(data = varImportanceMDG, aes(x = reorder(feature, MeanDecreaseGini), 
                                      y = MeanDecreaseGini,
                                      fill = MeanDecreaseGini)) + 
    geom_bar(stat = 'identity', color = "black", size = 0.1, show.legend = FALSE) + 
    labs(x = "Variable", y = "Importance (MDG)") +
    coord_flip() + 
    theme_bw() + 
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9,"OrRd")),
                         values = rescale(varImportanceMDG$MeanDecreaseGini, 
                                          to = c(0, 1))) + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) + 
    ggtitle("Variable Importance: Mean Decrease Gini")
  
  # save variable importance bar plot to image file 
  ggsave(filename = paste0(out_dir, outDescription,"varImpPlot_MDG_",
                           shapefile_layer_names$description[i],
                           ".png"))
  
  # make varImpPlot using the default dot plot in randomForest
  #randomForest::varImpPlot(rf_model,
  #                         main = shapefile_layer_names$description[i])
  
  
  # Ordered from highest MDGini to lowest
  write("\nVariable Importance, ranked by Mean Decrease Gini: \n", 
        rf_output_file, append=TRUE)
  capture.output(varImportanceMDG,
                 file = rf_output_file,
                 append=TRUE)
  
  write("\n", rf_output_file, append=TRUE)
  
  # RECORD TOP N MOST IMPORTANT VARIABLES BASED ON MEAN DECREASE GINI 
  rfVarImp[i,(2+nVar):((1+nVar)+nVar)] <- varImportanceMDG$feature[1:nVar]
  
  # variable importance, ordered from highest MDA to lowest
  write("\nVariable Importance, ranked by MDA: \n", 
        rf_output_file, append=TRUE)
  capture.output(varImportanceMDA,
                 file = rf_output_file,
                 append=TRUE)
  
  print(paste0("Top ", as.character(nVar)," most important variables ranked by MDA"))
  print(varImportanceMDA$feature[1:nVar]) 
  # RECORD TOP N MOST IMPORTANT VARIABLES BASED ON MEAN DECREASE ACCURACY 
  rfVarImp[i,2:(nVar+1)] <- varImportanceMDA$feature[1:nVar]
  
  
  # # TO DO: keep the n most important variables and run the classification again?
  # print("KEEPING TOP 5 VARIABLES AND RUNNING RF AGAIN ")
  # if (keepMostImpVar == TRUE) {
  #   mostImpVar <- rownames(varImp[1:5,])
  #   features2 <- features %>% dplyr::select(taxonID, c(mostImpVar))
  #   # run RF model again, this time with reduced features
  #   set.seed(104)
  #   rf_model2 <- randomForest::randomForest(as.factor(features2$taxonID) ~ .,
  #                                          data=features2,
  #                                          importance=TRUE,
  #                                          ntree=ntree) # ntree is number of trees to grow
  #   print(rf_model2)
  # }
  # initial test shows that this does not significantly improve the accuracy
  
  
  # INTERSPECIES VARIABLE COMPARISON BOXPLOTS -------------------------------
  
  if (createBoxplots == TRUE){
    
    # select the top n most important variables based on MDA for the current model 
    vars <- as.character(rfVarImp[i,2:(nVar+1) ])
    # subset the features data frame to obtain just the columns needed 
    test_features <- features %>% dplyr::select(c(taxonID, vars))
    # "gather" the wide data to make it longer.
    # each row represents a value. There is a separate column to indicate
    # species and also which feature (lidar, hs, rgb, etc...) for grouping in the plot
    features_gathered <- tidyr::gather(test_features, "feature", "value",-taxonID)
    # make a new column where each feature name is a factor, so the plots
    # can be arranged in order of most important to least important 
    features_gathered$feature_ordered = factor(features_gathered$feature, levels=vars)
    
    # create the multi-plot
    g <- ggplot(data = features_gathered, aes(x = taxonID, y = value)) + 
      geom_boxplot(aes(fill = taxonID),outlier.size = 0.2) + 
      facet_wrap(. ~ feature_ordered,scales='free',ncol=3) + 
      scale_fill_brewer(palette = "Spectral") + 
      theme_bw() + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.title = element_text(size = 10)) +
      ggtitle(paste0("Interspecies boxplot comparison of MDA most important features: \n",
                     shapefile_layer_names$description[i])) 
    ggsave(g, filename = paste0(out_dir,outDescription,"boxplot_impVarMDA_",shapefile_layer_names$description[i],".png"),
           width = 6, height = 5, dpi = 300, units = "in", device='png')
    
    
    # MDG
    # select the top n most important variables based on MDA for the current model 
    vars <- as.character(rfVarImp[i,(2+nVar):(1+2*nVar)])
    # subset the features data frame to obtain just the columns needed 
    test_features <- features %>% dplyr::select(c(taxonID, vars))
    # "gather" the wide data to make it longer.
    # each row represents a value. There is a separate column to indicate
    # species and also which feature (lidar, hs, rgb, etc...) for grouping in the plot
    features_gathered <- tidyr::gather(test_features, "feature", "value",-taxonID)
    # make a new column where each feature name is a factor, so the plots
    # can be arranged in order of most important to least important 
    features_gathered$feature_ordered = factor(features_gathered$feature, levels=vars)
    
    # create the multi-plot
    g2 <- ggplot(data = features_gathered, aes(x = taxonID, y = value)) + 
      geom_boxplot(aes(fill = taxonID),outlier.size = 0.2) + 
      facet_wrap(. ~ feature_ordered,scales='free',ncol=3) + 
      scale_fill_brewer(palette = "Spectral") + 
      theme_bw() + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.title = element_text(size = 10)) +
      ggtitle(paste0("Interspecies boxplot comparison of MDG most important features: \n",
                     shapefile_layer_names$description[i]))
    
    ggsave(g2, filename = paste0(out_dir,outDescription,"boxplot_impVarMDGini_",shapefile_layer_names$description[i],".png"),
           width = 6, height = 5, dpi = 300, units = "in", device='png')
    
  }
  
  
  write("\n\n------------------------------\n\n", rf_output_file, append=TRUE)
  
}

# close the text file
close(rf_output_file)

# write the accuracy summary data frame to file 
write.csv(rfAccuracies, paste0(out_dir, outDescription, "rfAccuraciesOOB.csv"))
# write the variable importance summary to file 
write.csv(rfVarImp, paste0(out_dir, outDescription, "rfVarImp.csv"))


end_time <- Sys.time()
print("Elapsed time: ")
print(end_time-start_time)




# ACCURACY TABLE  ---------------------------------------------------------

# create a nice table to summarize the model accuracies 
#https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html 
library(kableExtra)
rfAccuracies %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover","condensed"))


# Overall Accuracy Bar Plots ----------------------------------------------

# Accuracies displayed as bar plot 
oa_oob <- c(70.1, 73.1, 71.2, 61.0, 68.3, 67.4)
oa_indval <- c(56.2, 68.5, 68.5, 54.3, 65.4, 63.0)
ref_labs <- c("all points", "all polygons \nhalf diameter", "all polygons \nmax diameter",
              "neon_veg points", "neon_veg polygons \nhalf diameter", "neon_veg polygons \nmax diameter")

# convert from wide to long
# for a grouped bar plot for both kinds of error 
# USING ALL OBSERVATIONS PER REFERENCE DATA SET
acc_df <- data.frame(ref_labs, OOB=oa_oob, IndVal=oa_indval)
acc_df_long <- tidyr::gather(acc_df, oa_type, oa, OOB:IndVal)

acc_df_long$oa_type<-as.factor(acc_df_long$oa_type)

ggplot(data=acc_df_long, aes(x=ref_labs, y=oa, fill=oa_type)) +
  geom_bar(stat="identity", color="black", size = 0.25,
           position=position_dodge())+
  theme_bw() + 
  coord_cartesian(ylim=c(50, 75)) + 
  labs(x = "Reference Data Set", y = "OA [%]\n", 
       fill = "OA metric")+ 
  #geom_text(stat='identity', aes(label=oa), vjust=-0.5) + 
  geom_text(aes(label=oa), position=position_dodge(width=0.9), vjust=-0.25)+
  ggtitle("Overall Accuracy (OA) Comparison") + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_brewer(palette = "Greys")

ggsave(filename = paste0(out_dir,"OverallAccuracyComparisonBarPlot.png"),
       width = 8, height = 5, dpi = 500, units = "in", device='png')


##############################################
# OA comparison when SAMPLE SIZE IS REDUCED 

oa_oob <- c(61.5, 67.9, 63.0, 61.0, 68.3, 67.4)
oa_indval <- c(56.2, 65.4, 65.4, 54.3, 65.4, 63.0)
ref_labs <- c("all points", "all polygons \nhalf diameter", "all polygons \nmax diameter",
              "neon_veg points", "neon_veg polygons \nhalf diameter", "neon_veg polygons \nmax diameter")

# convert from wide to long
# for a grouped bar plot for both kinds of error 
# USING ALL OBSERVATIONS PER REFERENCE DATA SET
acc_df <- data.frame(ref_labs, OOB=oa_oob, IndVal=oa_indval)
acc_df_long <- tidyr::gather(acc_df, oa_type, oa, OOB:IndVal)

acc_df_long$oa_type<-as.factor(acc_df_long$oa_type)

ggplot(data=acc_df_long, aes(x=ref_labs, y=oa, fill=oa_type)) +
  geom_bar(stat="identity", color="black", size = 0.25,
           position=position_dodge())+
  theme_bw() + 
  coord_cartesian(ylim=c(50, 75)) + 
  labs(x = "Reference Data Set", y = "OA [%]\n", 
       fill = "OA metric")+ 
  #geom_text(stat='identity', aes(label=oa), vjust=-0.5) + 
  geom_text(aes(label=oa), position=position_dodge(width=0.9), vjust=-0.25)+
  ggtitle("Overall Accuracy (OA) Comparison with reduced sample size") + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_brewer(palette = "Greys")

ggsave(filename = paste0(out_dir,"OverallAccuracyComparisonBarPlot-ReducedSampleSize.png"),
       width = 8, height = 5, dpi = 500, units = "in", device='png')



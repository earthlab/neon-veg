# This script downloads NEON AOP remote sensing data 
# that will become descriptive features for species classification.
#
# Each data product is downloaded to a folder named by its 
# data product ID. This script also moves the downloaded files into
# data_raw/SITE_YEAR/ to more easily find the data later. 


# Download AOP data -------------------------------------------------------

# Download NEON AOP remote sensing tiles that cover the same area as the 
# in-situ woody vegetation data. We can specify specific coordinates using
# the tree locations from the woody vegetation structure data. 


# CHM ---------------------------------------------------------------------

dp_chm <- "DP3.30015.001"

neonUtilities::byTileAOP(
  dpID = dp_chm
  ,site = site_code
  ,year = aop_data_year
  ,savepath = dir_data_raw
  ,easting = veg_coordinates$eastings
  ,northing = veg_coordinates$northings
  ,buffer = buffer_val)

# move CHM files
move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_chm
                     ,dp_name = "chm", file_pattern = "*CHM.tif$")


# aspect_slope ------------------------------------------------------------

dp_aspect_slope <- "DP3.30025.001"

neonUtilities::byTileAOP(
  dpID = dp_aspect_slope
  ,site = site_code
  ,year = aop_data_year
  ,savepath = dir_data_raw
  ,easting = veg_coordinates$eastings
  ,northing = veg_coordinates$northings
  ,buffer = buffer_val)

# move aspect files
move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_aspect_slope
                      ,dp_name = "aspect", file_pattern = "*aspect.tif$"
                      ,delete_orig = FALSE)

# move slope files; delete originally downloaded aspect and slope files 
move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_aspect_slope
                      ,dp_name = "slope", file_pattern = "*slope.tif$"
                      ,delete_orig = TRUE)


# RGB ---------------------------------------------------------------------

dp_rgb <- "DP3.30010.001"

neonUtilities::byTileAOP(
  dpID = dp_rgb
  ,site = site_code
  ,year = aop_data_year
  ,savepath = dir_data_raw
  ,easting = veg_coordinates$eastings
  ,northing = veg_coordinates$northings
  ,buffer = buffer_val)

# move rgb files
move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_rgb
                      ,dp_name = "rgb", file_pattern = "*image.tif$"
                      ,delete_orig = TRUE)


# Vegetation indices ------------------------------------------------------

dp_veg_indices <- "DP3.30026.001"

neonUtilities::byTileAOP(
  dpID = dp_veg_indices
  ,site = site_code
  ,year = aop_data_year
  ,savepath = dir_data_raw
  ,easting = veg_coordinates$eastings
  ,northing = veg_coordinates$northings
  ,buffer = buffer_val)

move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_veg_indices
                      ,dp_name = "veg_indices"
                      ,file_pattern = "*VegIndices.zip$"
                      ,delete_orig = TRUE, unzip = TRUE)



# Hyperspectral reflectance -----------------------------------------------

dp_hs_refl <- "DP3.30006.001"

neonUtilities::byTileAOP(
  dpID = dp_hs_refl
  ,site = site_code
  ,year = aop_data_year
  ,savepath = dir_data_raw
  ,easting = veg_coordinates$eastings
  ,northing = veg_coordinates$northings
  ,buffer = buffer_val)


move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_hs_refl
                      ,dp_name = "hyperspectral"
                      ,file_pattern = "*reflectance.h5$"
                      ,delete_orig = TRUE)






# Download all data for a given site, year, data product ------------------

# To download all AOP data for specified site, year, product 
#neonUtilities::byFileAOP(
#  dpID = dp_chm
#  ,site = site_code
#  ,year = aop_data_year
#  ,check.size = TRUE
#  ,savepath = dir_data_raw)




# Looping through data product downloads ----------------------------------

# List of Data Product IDs to download
# "DP3.30015.001" - "chm"               LiDAR-derived Canopy Height Model
# "DP3.30025.001" - "aspect_slope"      LiDAR-derived Aspect and Slope
# "DP3.30010.001" - "rgb"               High-resolution multispectral camera RGB 
# "DP3.30026.001" - "veg_indices"       Spectrometer-derived vegetation indices 
# "DP3.30006.001" - "hs_refl"           Spectometer-derived surface reflectance

# data_products <- data.frame(id = c("DP3.30015.001"
#                                    ,"DP3.30025.001"
#                                    ,"DP3.30010.001"
#                                    ,"DP3.30026.001"
#                                    ,"DP3.30006.001")
#                             ,name = c("chm"
#                                       ,"aspect_slope"
#                                       ,"rgb"
#                                       ,"veg_indices"
#                                       ,"hs_refl"))


# Can this be turned into a loop?
# How to make R wait for user input (y/n) and wait for the files to download? 

# for (i in nrow(data_products)){
#   
#   dp_id <- as.character(data_products$id[i])
#   dp_name <- as.character(data_products$name[i])
#   
#   print(paste("Downloading data product ID", dp_id))
#   
#   neonUtilities::byTileAOP(
#     dp_id = dp_id
#     ,site = site_code
#     ,year = aop_data_year
#     ,savepath = dir_data_raw
#     ,easting = veg_coordinates$eastings
#     ,northing = veg_coordinates$northings
#     ,buffer = buffer_val)
  
} 
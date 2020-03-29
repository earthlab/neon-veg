# This script is where packages are installed and parameters are defined for 
# the other scripts within this workflow. 

# About neonUtilities: https://www.neonscience.org/neonDataStackR
# About geoNEON: https://github.com/NEONScience/NEON-geolocation/tree/master/geoNEON

# Load necessary packages. 
# If R says "Error in loadNamespace(name) : there is no package called 'packageName',
# Use install.packages("package_name") unless otherwise specified. 
library(neonUtilities) 
library(geoNEON)
#install.packages("devtools")
#library(devtools)
#install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
library(dplyr)
library(stringr)
library(sf)
# following packages needed for the clip_overlap custom function
library(sp) 
library(raster)
library(rgeos)
library(tools)
library(randomForest) 
library(gridExtra)


# Load custom supporting functions 
source("00-supporting_functions.R")

# Create folders to store the raw and processed data
check_create_dir("data")
dir_data_raw <- file.path("data","data_raw")
dir_data_out <- file.path("data","data_output")
check_create_dir(dir_data_raw)
check_create_dir(dir_data_out)
# Create a folder to store outputs for analysis 
dir_results <- file.path("results")
check_create_dir(dir_results)



# Download in-situ woody vegetation data ----------------------------------

# USER-DEFINED-INPUTS
site_code <- "NIWO"         # Four-digit NEON site code, character string type.
veg_start_date <- "2016-01" # Starting date and ending dates of woody veg data 
veg_end_date <- NA          # to download, using the "YYYY-MM" format. 
                            # Set to NA for earliest/latest data possible.

# Define parameters for the in-situ woody plant vegetation data to download.
# Let's load the in-situ Woody Vegetation Structure data straight into R.
# Specify the NEON site and starting/ending date(s) for the data. 
# A message in the console will display the total file size to be downloaded.
# Proceed by typing "y" and pressing Enter. 
# The woody vegetation structure tables are 
veg_raw <- neonUtilities::loadByProduct(dpID = "DP1.10098.001"   
                                        ,site = site_code              
                                        ,startdate = veg_start_date      
                                        ,enddate = veg_end_date
                                        ,package = "basic"          
                                        ,check.size = T)


# Create features (points or polygons) for each tree ----------------------

# Points, polygons with the maximum crown diameter, and polygons with 
# half the maximum crown diameter per tree are created and written
# to .shp files in the data/data_output directory. 

source("02-create_tree_features.R")


# Process tree features ---------------------------------------------------

# USER-DEFINED INPUT
px_thresh <- 2 # area threshold, integer count of of 1m by 1m pixels, 
               # to filter out small tree crowns in the processing workflow

source("03-process_tree_features.R")


# Download airborne remote sensing data -----------------------------------

# Define parameters for the AOP remote sensing data to download
aop_data_year <- "2017"   # four-digit year in character string "YYYY" for AOP imagery
buffer_val <- 5 #[m]       # integer buffer size around tree coordinates

source("04-download_aop_imagery.R")


# Prepare AOP imagery for training set creation ---------------------------

# This script combines all airborne remote sensing data layers per square-km 
# tile into a single multi-layer raster, saved as an .rds file. 
source("05-prep_aop_imagery.R")

# read the text file with wavelength data for subsequent steps
wavelengths = as.numeric(unlist(read.table(file.path(dir_data_out,
                                                     "wavelengths.txt"),
                                           skip = 1,
                                           col.names = 'wavelength')))


# Plot AOP imagery --------------------------------------------------------

# USER-DEFINED-INPUT
# Specify which tile to plot. 452000_4432000 has a visible road at NIWO site. 
east_north_string <- "452000_4432000"

source("06-plot_aop_imagery.R")


# Create training data for species classification -------------------------

# USER-DEFINED INPUT 
# shapefile containing tree points or polygons

shapefile_filename <- file.path(dir_data_out, "veg_points_w_height_diam.shp")
source("07-extract_training_features.R")
shapefile_filename <- file.path(dir_data_out, "veg_polygons_half_diam.shp")
source("07-extract_training_features.R")
shapefile_filename <- file.path(dir_data_out, "veg_polygons_max_diam.shp")
source("07-extract_training_features.R")
shapefile_filename <- file.path(dir_data_out, "veg_polys_half_diam_clipped_overlap.shp")
source("07-extract_training_features.R")
shapefile_filename <- file.path(dir_data_out, "veg_polys_max_diam_clipped_overlap.shp")
source("07-extract_training_features.R")
shapefile_filename <- file.path(dir_data_out, "veg_points_half_diam_clipped_overlap.shp")
source("07-extract_training_features.R")



# Classify species --------------------------------------------------------
# Train random forest model and assess species classification accuracy.
# Model outputs will be written to a folder within data/data_output 
# starting with "rf_" followed by a description of each shapfile 
# containing points or polygons per tree. 

# FOR THE ANALYSIS, extract the validation set of pixels -----------
independentValidationSet <- TRUE # if TRUE, keep separate set for validation

# randomly select this amount of data for training, use the rest for validation
percentTrain <- 0.8 

if(independentValidationSet){
  # randomly select <percentTrain> of data from the clipped half diameter 
  # polygons for training, use the remaining samples for validation.
  # keep track of which pixelNumber, easting, northing
  # that these pixels correspond to. 
  
  # read the file with spectra to sample randomly for validation trees.
  # Clipped half diameter polygons, in this case. 
  validationSourceFilename <- file.path(dir_data_out, 
          "veg_polys_half_diam_clipped_overlap-extracted_features.csv")
  
  # remove any spectra with a height of 0 and remove any factors
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
  
  # in the subsequent species classification step, when each set of 
  # training features is prepared for the classifier, any spectra with
  # pixelNumer_eastingIDs_northingIDs that belong in the validation set, 
  # "valIDs" will be removed. 
  
} 


# FOR THE ANALYSIS, determine the individual IDs of each clipped set.  
# To remove the sample size bias, if TRUE this filters down each of the 
# raw NEON shapefile data sets to only contain the individualIDs 
# present in the clipped set. This actually only makes sense for the polygon sets.
neonvegIDsForBothShapefiles <- TRUE





# Train random forest classifier for each of the training sets: 

# all tree point locations with height and crown diameter measurements 
shapefile_filename <- file.path(dir_data_out, 
                                "veg_points_w_height_diam.shp")
source("08-classify_species.R")

# all tree polygons, circular with half the maximum crown diameter per tree
shapefile_filename <- file.path(dir_data_out, 
                                "veg_polygons_half_diam.shp")
source("08-classify_species.R")

# all tree polygons, circular with the maximum crown diameter per tree
shapefile_filename <- file.path(dir_data_out, 
                                "veg_polygons_max_diam.shp")
source("08-classify_species.R")

# tree polygons, circular with half the maximum crown diameter per tree
# after being clipped and thresholded using experimental processing workflow
shapefile_filename <- file.path(dir_data_out, 
                                "veg_polys_half_diam_clipped_overlap.shp")
source("08-classify_species.R")

# tree polygons, circular with the maximum crown diameter per tree
# after being clipped and thresholded using experimental processing workflow
shapefile_filename <- file.path(dir_data_out, 
                                "veg_polys_max_diam_clipped_overlap.shp")
source("08-classify_species.R")

# tree points, corresponding to polygons with half the maximum crown diameter per tree
# after being clipped and thresholded using experimental processing workflow
shapefile_filename <- file.path(dir_data_out, 
                                "veg_points_half_diam_clipped_overlap.shp")
source("08-classify_species.R")

# Accuracy assessment / comparison ----------------------------------------

# list the shapfile names corresponding to directories with models to assess
dirs_to_assess <- c("veg_points_w_height_diam.shp"
                    ,"veg_polygons_half_diam.shp"
                    ,"veg_polygons_max_diam.shp"
                    ,"veg_points_half_diam_clipped_overlap.shp"
                    ,"veg_polys_half_diam_clipped_overlap.shp"
                    ,"veg_polys_max_diam_clipped_overlap.shp")

source("09-assess_accuracy.R")



# Plots ----------------------------------------------------------------------


# Create study area map - highlight the NIWO location within Colorado, USA
# Create ribbon plots - display spectral signatures extracted from AOP data per species 
source("10-create_figures.R")


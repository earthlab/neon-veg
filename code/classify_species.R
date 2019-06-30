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
source("supporting_functions.R")

# code for NEON site 
site_code <- "NIWO"

# define the "bad bands" wavelength ranges in nanometers, where atmospheric 
# absorption creates unreliable reflectance values. 
bad_band_window_1 <- c(1340, 1445)
bad_band_window_2 <- c(1790, 1955)

# taxon ID's to predict
taxonList <- c("ABLAL","PICOL","PIEN","PIFL2")

# define the output directory. If it doesn't exist already, create it.
check_create_dir('../output/') # create top level "output" directory
out_dir <- paste0('../output/', site_code, '/')
check_create_dir(out_dir) # create output folder for site
check_create_dir(paste0(out_dir,"figures/")) # create output figures folder

# shapefile set(s) to test --------------------------------------------------

# shapefiles were generated using the neon_veg workflow
# where each entry (point or polygon) corresponds to a tree in the woody 
# vegetationin-situ field NEON data set.

# directory with shapefiles to test (tree stem locations and crown polygons)
shapefile_dir <- paste0("output/", site_code, "/")

# define the subdirectory (destination or dsn when reading shapefile)
# and the layer name (the filename before the .shp extension) 
# in a vector for each of the shapefile scenarios. 

# all stem points at NIWO that also have crown diameter measurements. 
# this includes multi-bole entries, which result in duplicated spectra. 
allStems_layer <- c("allStems",
                    "shapefiles_maxDiameter/",
                    "mapped_stems_with_crown_diameter")

# polygons with max crown diameter, one generated for each stem point at NIWO.
# just as for the "allStems" layer above, the shapfile used here
# includes multi-bole entries, since this is how the data come when
# downloaded straight from the data portal. 
allPolygons_maxDiameter_layer <- c("allPolygons_maxDiameter",
                                   "shapefiles_maxDiameter/",
                                   "polygons_all")

# polygons with half max diameter, one generated for each stem point at NIWO. 
allPolygons_halfDiameter_layer <- c("allPolygons_halfDiameter",
                                    "shapefiles_halfDiameter/",
                                    "polygons_all")

# neon_veg workflow stems, corresponding to polygons generated with max diameter
neonvegStems_maxDiameter_layer <- c("neonvegStems_maxDiameter",
                                    "shapefiles_maxDiameter/",
                                    "mapped_stems_final")

# neon_veg workflow polygons generated with half the max diameter  
neonvegPolygons_halfDiameter_layer <- c("neonvegPolygons_halfDiameter",
                                        "shapefiles_halfDiameter/",
                                        "polygons_clipped_overlap")

# neon_veg workflow polygons generated with max diameter 
neonvegPolygons_maxDiameter_layer <- c("neonvegPolygons_maxDiameter",
                                       "shapefiles_maxDiameter/",
                                       "polygons_clipped_overlap")

# create a data frame containing the description, directory, and shapefile name
# for each of the shapefile scenarios to be tested 
shapefileLayerNames <- as.data.frame(rbind(allStems_layer,
                                           allPolygons_halfDiameter_layer,
                                           allPolygons_maxDiameter_layer,
                                           neonvegStems_maxDiameter_layer,
                                           neonvegPolygons_halfDiameter_layer,
                                           neonvegPolygons_maxDiameter_layer),
                                     stringsAsFactors = FALSE) %>% 
  `colnames<-`(c("description", "dsn", "layer")) 
rownames(shapefileLayerNames) <- 1:nrow(shapefileLayerNames)

print("shapefile layers to test:")
shapefileLayerNames

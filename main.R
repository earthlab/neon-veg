# load packages 
library(devtools)
library(geoNEON)
library(geosphere)

# load functions written in external files. 
# setwd("~/github/neon-veg")
source("locate_woody_veg.R")
source("woody_df_to_shp.R")
source("merge_vst_tables.R")


# define path to NEON l1 Woody Vegetation data 
woody_path <- "~/Documents/earth-analytics/data/SJER_2017/NEON_struct-woody-plant/NEON.D17.SJER.DP1.10098.001.2015-04.basic.20171016T160520Z/"  

# mapping and tagging table (contains stem locations)
woody_mapping_path <- paste(woody_path, 
                       list.files(path = woody_path, pattern = "mappingandtagging"), 
                       sep="")

# apparent individual table (contains crown shape)
woody_individual_path <- paste(woody_path, 
                          list.files(path = woody_path, pattern = "apparentindividual"), 
                          sep="")


# load in situ stem locations table
woody <- read.csv(woody_mapping_path)

# calculate mapped UTM locations of plants from distance/azimuth
woody_utm <- locate_woody_veg(woody)

# load "vst_apparentindividual" table
woody_individual <- read.csv(woody_individual_path)

# match mapped stems from "vst_mappingandtagging" with structure data from "vst_apparentindividual"
woody_vst <- merge_vst_tables(woody_utm, woody_individual)

# create circular polygon for each stem based on maxCanopyDiameter
woody_polygons <- woody_df_to_shp()
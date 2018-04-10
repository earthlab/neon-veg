# load packages 
library(devtools)
library(geoNEON)
library(sp)           # CRS
library(swfscMisc)    # circle.polygon
library(rgdal)      # writeOGR

# load functions written in external files. 
# setwd("~/github/neon-veg")
source("locate_woody_veg.R")
source("woody_df_to_shp.R")
source("merge_vst_tables.R")
source("get_vst_crs.R")

# define path to NEON l1 Woody Vegetation data 
main_path <- "~/Documents/earth-analytics/data/SJER_2017/NEON_struct-woody-plant/"


# loop through folders of field data with different dates
dirs <- list.dirs(path = main_path )
dirs <-dirs[ grepl("NEON.D17.SJER", dirs) ]

first_loop <- 1 # loop counter
for (woody_path in dirs){
  
  # mapping and tagging table (contains stem locations)
  woody_mapping_path <- paste(woody_path, 
                              list.files(path = woody_path, 
                                         pattern = "mappingandtagging"), 
                              sep="/")
  
  # apparent individual table (contains height and crown diameter)
  woody_individual_path <- paste(woody_path, 
                                 list.files(path = woody_path, 
                                            pattern = "apparentindividual"), 
                                 sep="/")
  
  # load situ stem locations table; 
  # calculate mapped UTM locations of plants from distance/azimuth
  woody <- read.csv(woody_mapping_path)
  woody_utm <- locate_woody_veg(woody)
  
  # load "vst_apparentindividual" table
  woody_individual <- read.csv(woody_individual_path)
  
  # match mapped stems from "vst_mappingandtagging" with structure data from "vst_apparentindividual"
  woody_vst <- merge_vst_tables(woody_utm, woody_individual)
  
  # combine woody veg struct data to a single data frame 
  if (first_loop == 1){
    woody_all <- woody_vst
    first_loop <- 0
    
  } else {
    woody_all <- rbind(woody_all, woody_vst)
  }
}

# create coordinate reference system object based on
# UTM zone info in the "vst_plotperyear" table
crs <- get_vst_crs(woody_path)

# create circular polygon for each stem based on maxCanopyDiameter
woody_polygons <- woody_df_to_shp(df=woody_all, 
                                  coord_ref=crs,
                                  shrink=1,
                                  num_sides = 8,
                                  shp_filename = "test_sjer_polygons")

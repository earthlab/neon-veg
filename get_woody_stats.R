# get woody veg data statistics
# such as, how many trees are multi-bole? 
# the apparent individual data table contains a column called growthForm


# setup -------------------------------------------------------------------

working.dir <- "~/github/neon-veg/"
setwd(working.dir)

main.path <- "SJER/woody_veg" 

# specify output directory path and filename of output shapefile to be written
out.dir <- "output/"

# load required packages
library(geoNEON)
source("locate_woody_veg.R")
source("merge_vst_tables.R")


# get stats ---------------------------------------------------------------

# get list of directories containing NEON data 
dirs <- list.dirs(path = main.path )
dirs <- dirs[grepl("NEON.D", dirs)]


first.loop <- 1 # loop counter
for (woody.path in dirs) {
  
  # mapping and tagging table (contains stem locations)
  woody.mapping.path <- paste(woody.path, 
                              list.files(path = woody.path, 
                                         pattern = "mappingandtagging"), 
                              sep = "/")
  
  # apparent individual table (contains height and crown diameter)
  woody.individual.path <- paste(woody.path, 
                                 list.files(path = woody.path, 
                                            pattern = "apparentindividual"), 
                                 sep = "/")
  
  # load situ stem locations table; 
  # calculate mapped UTM locations of plants from distance/azimuth
  woody <- read.csv(woody.mapping.path)
  woody.utm <- locate_woody_veg(woody)
  
  # load "vst_apparentindividual" table
  woody.individual <- read.csv(woody.individual.path)
  
  
  # merge mapped stems from "vst_mappingandtagging" (x) with structure data 
  # from "vst_apparentindividual" (y) based on individualID 
  woody.vst <- merge(x = woody.utm, 
                     y = woody.individual,
                     by = "individualID")
  # create a new column for the "date" of collection from apparent_individual
  # keep most recent entry for multiple entries with same individual ID
  woody.vst$date <- woody.vst$date.y 
  woody.vst <- woody.vst %>% select(-c(date.x, date.y))
                  group_by(individualID) %>%
                  slice(which.max(as.Date(date)))
                  
  # count the number of "complete entries" with crown diameter and height 
  complete.cases <- woody.vst[complete.cases(woody.vst$maxCrownDiameter) & 
                              complete.cases(woody.vst$height),] %>% 
    select(individualID, height, maxCrownDiameter)
  
  print(paste(as.character(nrow(complete.cases)), 
              "trees have height and crown diameter measurements."))
           
                  
  # how many trees have crown diameter but NO height? 
  height.cases <- woody.vst[complete.cases(woody.vst$height),]
  
  print(paste(as.character(nrow(height.cases)), 
              "trees have height measurements."))
  
  
  # how many trees have height but NO crown diameter? 
  crown.diam.cases <- woody.vst[complete.cases(woody.vst$maxCrownDiameter),]
  
  print(paste(as.character(nrow(crown.diam.cases)), 
              "trees have crown diameter measurements."))
  
  
  # how many trees have multi-bole growth forms? 
  multi.bole.cases <- 
  
  
  
  
  # the apparent individual data table contains a column called growthForm
  
  
  # combine woody veg structure data to a single data frame 
  if (first.loop == 1) {
    woody_all <- woody.vst
    woody_mapping_all <- woody.utm
    woody.individual_all <- woody.individual
    first.loop <- 0
    
  } else {
    woody_all <- rbind(woody_all, woody.vst)
    woody_mapping_all <- rbind(woody_mapping_all, woody.utm)
    woody.individual_all <- rbind(woody.individual_all, woody.individual)
  }
}





# this code is based off of the following NEON tutorials:
# Working with Hyperspectral Data in HDf5 Format in R
# http://neondataskills.org/HDF5/Imaging-Spectroscopy-HDF5-In-R/ 

# use the code below to install the rhdf5 library if it's not already installed.
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

# r Load `raster` and `rhdf5` packages and read hyperspectral data into R
library(rhdf5)
library(rgdal)
library(raster)
library(ggplot2)
library(tidyr)
library(sf)
library(dplyr)
library(data.table) 

# set working directory
setwd("~/github/neon-veg/")

# load custom functions
source("supporting_functions.R") #check_create_dir


# directory with hyperspectral .h5 files
h5.dir <- '~/github/neon-veg/data/NIWO/hyperspectral_reflectance/'

# read the polygon/points data 
site <- "NIWO"
polygon.path <- '~/github/neon-veg/data/NIWO/veg_polygons/'

# specify the output directory to write spectra to file 
out.dir.spectra <- '~/github/neon-veg/data/NIWO/output_spectra/'
check_create_dir(out.dir.spectra)

# read polygon file 
polygons <- rgdal::readOGR(dsn = polygon.path,
                           layer = "polygons_checked_overlap")

# read the tree locations 
tree.points <- rgdal::readOGR(dsn = polygon.path,
                              layer = "mapped_stems_final")

# convert polygons and tree locations to SF objects
polygons.sf <- sf::st_as_sf(polygons)
tree.points.sf <- sf::st_as_sf(tree.points) 

# isolate the tree location coordinates 
tree.coords <- tree.points.sf %>% 
  sf::st_coordinates() %>% 
  as.data.frame()

# add new columns for the tree location coordinates 
tree.points.sf$X <- tree.coords$X
tree.points.sf$Y <- tree.coords$Y

# add empty columns for the min and max coordinates for each polygon
polygons.sf$xmin <- NA 
polygons.sf$xmax <- NA 
polygons.sf$ymin <- NA 
polygons.sf$ymax <- NA 

# add the min, max X and Y values to each polygon for filtering 
for (i in 1:nrow(polygons.sf)) {
  polygons.sf$xmin[i] <- as.numeric(sf::st_bbox(polygons.sf$geometry[i])[1])
  polygons.sf$ymin[i] <- as.numeric(sf::st_bbox(polygons.sf$geometry[i])[2])
  polygons.sf$xmax[i] <- as.numeric(sf::st_bbox(polygons.sf$geometry[i])[3])
  polygons.sf$ymax[i] <- as.numeric(sf::st_bbox(polygons.sf$geometry[i])[4])
}

# merge the polygons with tree locations;
# rename the geometry columns to be more descriptive 
polygons.points <- merge(as.data.frame(polygons.sf),
                         as.data.frame(tree.points.sf)
                         [,c("indvdID", "X","Y","geometry")],
                         by="indvdID") %>% 
  dplyr::rename(geometry.polygon = geometry.x, 
                geometry.point = geometry.y)










# specify red, green, and blue wavelengths for RGB composite
rgb.bands <- c(669, 549, 474) 

# specify "bad band" windows, in nanometers
bad.band.window.1 <- c(1340, 1445)
bad.band.window.2 <- c(1790, 1955)


# loop through h5 files ---------------------------------------------------

# get the names of all h5 files to iterate through
h5.list <- list.files(path = h5.dir, full.names = TRUE)
h5.list <- h5.list[grepl("*.h5", h5.list)]

# loop through h5 files 
for (h5 in h5.list) {
  
  print(h5)
  
  # list the contents of the HDF5 file
  h5.struct <- rhdf5::h5ls(h5, all=T)
  
  # construct the string using "/Reflectance/Metadata/Coordinate_System",
  # without explicitly using a site code 
  crs.tag <- h5.struct$group[grepl("/Reflectance/Metadata/Coordinate_System", 
                                   h5.struct$group)][1] 
  
  # read coordinate reference system data
  crs.info <- rhdf5::h5read(h5, crs.tag)
  
  # convert "UTM" to lowercase "utm" for proper usage later
  crs.info$Proj4 <- CRS(chartr("UTM", "utm", crs.info$Proj4))
  
  # print the coordinate reference data
  crs.info
  
  # get attributes for the Reflectance dataset.
  # construct the string using /Reflectance/Reflectance_Data 
  refl.tag <- paste0(h5.struct$group[grepl("/Reflectance", 
                                           h5.struct$group)][1],
                     "/Reflectance_Data")
  
  refl.info <- rhdf5::h5readAttributes(h5,refl.tag)
  
  n.rows <- refl.info$Dimensions[1]
  n.cols <- refl.info$Dimensions[2]
  n.bands <- refl.info$Dimensions[3]
  
  # print dimensions 
  print(paste0("# Rows: ", as.character(n.rows)))
  print(paste0("# Columns: ", as.character(n.cols)))
  print(paste0("# Bands: ", as.character(n.bands)))
  
  # read the wavelengths of the hyperspectral image bands
  wavelength.tag <- paste0(h5.struct$group[grepl("/Reflectance/Metadata/Spectral_Data", 
                                                 h5.struct$group)][1],
                           "/Wavelength")
  wavelengths <- rhdf5::h5read(h5,
                               wavelength.tag)
  
  # define spatial extent: extract resolution and origin coordinates
  map.info <- unlist(strsplit(crs.info$Map_Info, 
                              split = ", "))
  res.x <- as.numeric(map.info[6])
  res.y <- as.numeric(map.info[7])
  x.min <- as.numeric(map.info[4])
  y.max <- as.numeric(map.info[5])
  
  # calculate the maximum X and minimum Y values 
  x.max <- (x.min + (n.cols*res.x))
  y.min <- (y.max - (n.rows*res.y))
  tile.extent <- raster::extent(x.min, x.max, y.min, y.max)
  
  print("tile extent: ")
  tile.extent
  
  # figure out which trees are within the current tile 
  polygons.in <- polygons.points %>% 
    dplyr::filter(xmin > tile.extent@xmin & 
                    xmax < tile.extent@xmax & 
                    ymin > tile.extent@ymin & 
                    ymax < tile.extent@ymax)
  
  print(paste0(as.character(nrow(polygons.in))," polygons in current tile"))
  
  # if no polygons are within the current tile, skip to the next one
  if (nrow(polygons.in)==0){
    print("no trees located within current tile")
    next
  }
  
  # read reflectance data for all bands
  refl <- rhdf5::h5read(h5,refl.tag,
                        index = list(1:n.bands, 1:n.cols, 1:n.rows))
  
  # view and apply scale factor to convert integer values to reflectance [0,1]
  # and data ignore value
  scale.factor <- refl.info$Scale_Factor
  data.ignore <- refl.info$Data_Ignore_Value
  refl[refl == data.ignore] <- NA 
  refl.scaled <- refl / scale.factor
  
  # ########## testing polygon extraction ###################
  # # create georeferenced raster using band 1 
  # refl.1 <- (refl.scaled[1,,]) # convert to matrix
  # refl.raster <- raster::raster(refl.1, crs = crs.info$Proj4)
  # extent(refl.raster) <- tile.extent
  # 
  # # extract reflectance values from raster pixels inside polygon boundary 
  # raster::extract(refl.raster, sf::as_Spatial(polygons.in$geometry.polygon[1]))
  # # finds the intersection between bounding boxes of two spatial objects 
  # test <- raster::intersect(refl.raster, sf::as_Spatial(polygons.in$geometry.polygon[5])); plot(test)
  # # finds the intersection between geometries
  # test <- rgeos::gIntersection(refl.raster, sf::as_Spatial(polygons.in$geometry.polygon[3])); plot(test)
  # #########################################################
  
  # plot RGB composite 
  # plot_hs_rgb(refl.scaled, 
  #             wavelengths,
  #             rgb.bands, 
  #             crs.info$Proj4, 
  #             tile.extent,
  #             plt=TRUE)
  
  # create data frame to store extracted spectra 
  spectra <- data.frame(matrix(NA, 
                               nrow = dim(wavelengths), 
                               ncol = nrow(polygons.in)))
  
  # extract hyperspectral reflectance at each tree location (point)
  for (i in 1:nrow(polygons.in)){
    # get the image coordinates of each tree location
    index.x <- round(polygons.in$X[i] - tile.extent@xmin)
    index.y <- round(polygons.in$Y[i] - tile.extent@ymin)
    # extract the spectra 
    spectra[,i] <- refl.scaled[,index.x, index.y]
  }
  
  # format the spectra for writing to file 
  spectra.write <- data.frame(individualID = polygons.in$indvdID,
                              scientificName = polygons.in$scntfcN,
                              taxonID = polygons.in$taxonID,
                              maxCrownDiameter = polygons.in$crownDm,
                              height = polygons.in$height,
                              X = polygons.in$X,
                              Y = polygons.in$Y,
                              t(spectra))
  
  # rename the columns with rounded wavelength values
  data.table::setnames(spectra.write, 
                       old = colnames(spectra.write[,-(1:7)]), 
                       new = as.character(round(wavelengths)))
  
  write.csv(spectra.write, file = paste0(out.dir.spectra,"spectral_reflectance_",
                                          as.character(tile.extent@xmin),"_",
                                          as.character(tile.extent@ymin),".csv"))  
}



# read and plot the spectra .csv files  -----------------------------------

# get a list of the .csv file per tile containing woody veg stemse
csvs <- list.files(out.dir.spectra, full.names = TRUE)
csvs <- csvs[grepl("*000.csv", csvs)]

# combine all .csv data into a single data frame 
for (i in 1:length(csvs)){
  print(csvs[i])
  csv <- read.csv(csvs[i])
  
  if(i==1){
    spectra.all <- csv
  } else {
    spectra.all <- rbind(spectra.all, csv)
  }
}

# remove the unneccessary column 
spectra.all <- spectra.all %>% dplyr::select(-X.1)

# write ALL the spectra to a single file 
write.csv(spectra.all,
          file=paste0(out.dir.spectra,site, "_spectral_reflectance_ALL.csv"))

# write the exact wavelengths to file for future use 
write.table(data.frame(wavelengths = wavelengths),
            paste(out.dir.spectra,"wavelengths.txt"),
            sep="\n",
            row.names=FALSE)

# read wavelengths if not previously created
#wavelengths = read.table(paste(out.dir.spectra,"wavelengths.txt"),
#                         sep="\n",
#                         skip = 1,
#                         col.names = 'wavelength')

# Plot spectra and color line by species ----------------------------------

# keep just the reflectance values along with individual ID and species
# for the first tree
spectra.1 <- spectra.all[1,1:ncol(spectra.all)] %>% 
                  dplyr::select(-c(X,Y))

# organize the data into a data frame where the first column contains
# the individual ID, second column contains scientific name, 
# third column contains wavelength, and fourth column contains
# reflectance per wavelength.  
spectra.plot <- tidyr::gather(spectra.1,
                              key = wavelength,
                              value = "reflectance",
                              -individualID,
                              -scientificName,
                              -taxonID,
                              -maxCrownDiameter,
                              -height)
# add the actual wavelenth values into the wavelengths column
spectra.plot$wavelength <- wavelengths


for (i in 2:nrow(spectra.all)){
  spectra.current <- spectra.all[i,1:ncol(spectra.all)] %>% 
                        dplyr::select(-c(X,Y))
  spectra.current <- tidyr::gather(spectra.current,
                                key = wavelength,
                                value = "reflectance",
                                -individualID,
                                -scientificName,
                                -taxonID,
                                -maxCrownDiameter,
                                -height)
  spectra.current$wavelength <- wavelengths
  
  spectra.plot <- rbind(spectra.plot, spectra.current)
  
}

# remove the first reflectance value 
spectra.plot <- spectra.plot[spectra.plot$wavelength > 385,]

# remove the bad bands 
remove.bands <- wavelengths[(wavelengths > bad.band.window.1[1] & 
                            wavelengths < bad.band.window.1[2]) | 
                           (wavelengths > bad.band.window.2[1] & 
                            wavelengths < bad.band.window.2[2])]

spectra.plot$reflectance[spectra.plot$wavelength %in% remove.bands] <- NA

# plot all spectra on single plot 
legend.labels = c('a','b','c','d')
ggplot(data = spectra.plot, aes(x = wavelength, y = reflectance, colour = scientificName)) +
  geom_line(aes(group = individualID)) + 
  labs(x = "wavelength (nm)", color = "species") + 
  ggtitle("Hyperspectral reflectance extracted per center pixel") 
          + scale_color_manual(labels = legend.labels, 
                     values = c("#d7191c", "#fdae61", "#abdda4", "#2b83ba")
  )


# Use facet_grid for a separate plot per species 
ggplot(data = spectra.plot, 
       aes(x = wavelength, 
           y = reflectance, 
           colour = taxonID)) +
  geom_line(aes(group = individualID), alpha = 0.5) + 
  facet_wrap(. ~ taxonID, ncol = 2) + 
  labs(x = "wavelength (nm)", color = "species") + 
  ggtitle("Hyperspectral reflectance extracted per center pixel") + 
  scale_color_manual(
                     values = c("#d7191c", "#fdae61", "#abdda4", "#2b83ba"))

# save plot to image file 
ggsave(paste0(out.dir.spectra,"spectral_reflectance_per_taxonID.png"))

# spectral distance calculations ------------------------------------------




# raster stack HS image ---------------------------------------------------

# create georeferenced raster using band 1 
r1 <- (refl.scaled[1,,]) # convert first band to matrix
r1 <- raster::raster(r1, crs = crs.info$Proj4)
extent(r1) <- tile.extent
# start the raster stack with first band 
s <- raster::stack(r1)

# loop through bands and create a giant rasterstack with 426 bands
for(b in 2:n.bands){
  print(b)
  
  # create raster with current band
  r <- (refl.scaled[b,,]) # convert to matrix
  r <- raster::raster(r, crs = crs.info$Proj4)
  extent(r) <- tile.extent
  
  # add additional band to the stack with the addLayer function
  s <- raster::addLayer(s, r)
  
}

# adjust the names for each layer in raster stack to correspond to wavelength
names(s) <- round(wavelengths)

# instead of using the sf version of the polygons, need to use SpatialPolygons 
# find which trees are within current tile 
tree.coords.in <- sp::point.in.polygon(tree.coords$X, tree.coords$Y,
  c(tile.extent@xmin, tile.extent@xmax, tile.extent@xmax, tile.extent@xmin),
  c(tile.extent@ymin, tile.extent@ymin, tile.extent@ymax, tile.extent@ymax))

# subset the SpatialPolygons for those within current tile
# ?????? Get single polygon from the SpatialPolygonsDataFrame
spatial.polygons.in <- polygons@polygons[tree.coords.in]
p <- spatial.polygons.in[1]

# extract the spectra within each polygon within current tile
extracted <- raster::extract(x = s,
                             y = p,
                             df = TRUE)









# extract spectra within each polygon -------------------------------------

# create a georeferenced raster using the first hyperspectral band (just for 
# the coordinates, to compare to each polygon)

band2Raster <- function(file, band, noDataValue, xMin, yMin, res, crs){
  #first read in the raster
  out<- h5read(f,"Reflectance",index=list(1:nCols,1:nRows,band))
  #Convert from array to matrix
  out <- (out[,,1])
  #transpose data to fix flipped row and column order 
  #depending upon how your data are formated you might not have to perform this
  #step.
  out <-t(out)
  #assign data ignore values to NA
  #note, you might chose to assign values of 15000 to NA
  out[out == myNoDataValue] <- NA
  
  #turn the out object into a raster
  outr <- raster(out,crs=myCrs)
  
  # define the extents for the raster
  #note that you need to multiple the size of the raster by the resolution 
  #(the size of each pixel) in order for this to work properly
  xMax <- xMin + (outr@ncols * res)
  yMin <- yMax - (outr@nrows * res)
  
  #create extents class
  rasExt  <- extent(xMin,xMax,yMin,yMax)
  
  #assign the extents to the raster
  extent(outr) <- rasExt
  
  #return the raster object
  return(outr)
}

# create a georeferenced raster using band 1 to find which pixels each polygon 
# intersects with

# get image indices of the pixels within each polygon

# slice the entire hyperspectral reflectance, then average the pixel values? 
  # save them? 

# OR, clip the entire hyperspectral raster using the extent of the tree, 
  #then convert the remaining (much smaller ) 





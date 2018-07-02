
# This code is based off of the following NEON tutorials:
# Working with Hyperspectral Data in HDf5 Format in R
# http://neondataskills.org/HDF5/Imaging-Spectroscopy-HDF5-In-R/ 

# use the code below to install the rhdf5 library if it's not already installed.
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

# r Load `raster` and `rhdf5` packages and read NIS data into R
library(rhdf5)
library(rgdal)
library(raster)
library(ggplot2)
library(tidyr)

# set working directory
setwd("~/github/neon-veg/")

# load custom functions
source("get_hs_band.R")
source("plot_hs_rgb.R")


# Define the file name to be opened
f <- '~/github/neon-veg/NIWO/hyperspectral_reflectance/NEON_D13_NIWO_DP3_451000_4432000_reflectance.h5'

# look at the HDF5 file structure 
h5.struct <- rhdf5::h5ls(f,
                  all=T)

# read coordinate reference data
crs.info <- rhdf5::h5read(f, 
                   "NIWO/Reflectance/Metadata/Coordinate_System")

# convert "UTM" to lowercase "utm" for proper usage later
crs.info$Proj4 <- chartr("UTM", "utm", crs.info$Proj4)

# print the coordinate reference data
crs.info

# get attributes for the Reflectance dataset
refl.info <- rhdf5::h5readAttributes(f,
                              "NIWO/Reflectance/Reflectance_Data")

# resolution --> the Dimension Labels are "Line, Sample, Wavelength"
# does this correspond to n.rows first?????? in the NEON tutorial,
# it was stated that the h5 reader reads cols first instead of rows.
# this won't be a problem for square rasters, but are the ones on the 
# edges not square? 
# Dimensions are displayed as n.bands x n.rows (?) x n.cols (?) in refl.info
n.rows <- refl.info$Dimensions[1]
n.cols <- refl.info$Dimensions[2]
n.bands <- refl.info$Dimensions[3]

# print dimensions 
print(paste0("# Rows: ", as.character(n.rows)))
print(paste0("# Columns: ", as.character(n.cols)))
print(paste0("# Bands: ", as.character(n.bands)))


# wavelengths 
wavelengths <- rhdf5::h5read(f,
                      "NIWO/Reflectance/Metadata/Spectral_Data/Wavelength")

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

# read reflectance data
# 1 band
band = 1 
refl <- rhdf5::h5read(f,"/NIWO/Reflectance/Reflectance_Data",
                 index = list(band, 1:n.cols, 1:n.rows))
# all bands
refl <- rhdf5::h5read(f,"/NIWO/Reflectance/Reflectance_Data",
                index = list(1:n.bands, 1:n.cols, 1:n.rows)) 

# view and apply scale factor to convert integer values to reflectance [0,1]
# and data ignore value
scale.factor <- refl.info$Scale_Factor
data.ignore <- refl.info$Data_Ignore_Value
refl[refl == data.ignore] <- NA 
refl.scaled <- refl / scale.factor


# plot single band based on wavelength  -------------------------------------------------------

# specify wavelength
wl <- 669
get_hs_band(refl.scaled, 
            wavelengths,
            wl,
            crs.info$Proj4,
            tile.extent,
            plt = TRUE)


# plot RGB composite ----------------------------------------------------
# specify indices of the chosen red, green, and blue bands to display 
rgb.bands <- c(669, 549, 474) 
plot_hs_rgb(refl.scaled, 
            wavelengths,
            rgb.bands, 
            crs.info$Proj4, 
            tile.extent,
            plt=TRUE)


# plot reflectance spectra ------------------------------------------------

# select x,y pixel  
x <- 300
y <- 860
b <- 1:n.bands

spectrum <- cbind.data.frame(wavelength = wavelengths, 
                             refl = refl.scaled[x,y,b])

ggplot(data=spectrum, aes(x=wavelength, y = refl)) + 
    geom_line() + 
  labs(x = "Wavelength (nm)", y = "Reflectance") + 
  ggtitle("Spectral reflectance")



# clip using tree stem location points ------------------------------------

# test for a single polygon; plot RGB 

# polygon shapefiles to read
polygon.path <- "NIWO/output/"

# read polygon file 
polygons <- rgdal::readOGR(dsn = polygon.path,
                           layer = "polygons_checked_overlap")

# read the tree locations 
tree.points <- rgdal::readOGR(dsn = polygon.path,
                                 layer = "mapped_stems_final")

# use a plot boundary instead of tile extent 
clip.extent <- round(extent(c(451365.27,  # xmin
              451405.39,  # xmax
              4432738.62, # ymin
              4432778.74))) # ymax

# convert polyogns to SF object
library(sf)
polygons.sf <- sf::st_as_sf(polygons)

# convert tree locations to SF object
tree.points.sf <- sf::st_as_sf(tree.points) 
tree.coords <- tree.points.sf %>% 
  sf::st_coordinates() %>% 
  as.data.frame()
tree.points.sf$X <- tree.coords$X
tree.points.sf$Y <- tree.coords$Y

# add empty columns for the min and max coordinates 
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
                            as.data.frame(tree.points.sf)[,c("indvdID", "X","Y","geometry")],
                            by="indvdID") %>% 
  dplyr::rename(geometry.polygon = geometry.x, 
                geometry.point = geometry.y)
           
# filter the polygons, keep only the ones inside the current tile 
polygons.in <- polygons.points %>% dplyr::filter(xmin > clip.extent@xmin & 
                                          xmax < clip.extent@xmax & 
                                          ymin > clip.extent@ymin & 
                                          ymax < clip.extent@ymax)

# converted PYTHON to R: https://www.neonscience.org/neon-aop-hdf5-py 
# We will load the function calc_clip_index, which reads in a dictionary of 
# the spatial extent of the clipped region of interest, and a dictionary of 
# the full extent of the array you are subsetting, and returns the pixel 
# indices corresponding to the full flightline array.

full.extent <- tile.extent

calc_clip_index <- function(clip.extent, full.extent, x.scale=1, y.scale=1){
  
  h5rows = full.extent@ymax - full.extent@ymin
  h5cols = full.extent@xmax - full.extent@xmin

  ind.extent = {}
  ind.extent['xMin'] = round((clip.extent@xmin - full.extent@xmin)/x.scale)
  ind.extent['xMax'] = round((clip.extent@xmax - full.extent@xmin)/x.scale)
  ind.extent['yMax'] = round(h5rows - (clip.extent@ymin - full.extent@ymin)/x.scale)
  ind.extent['yMin'] = round(h5rows - (clip.extent@ymax - full.extent@ymin)/y.scale)

return(ind.extent)

}

subset.index = calc_clip_index(clip.extent,
                               tile.extent)
subset.index

# subset the reflectance array 
refl.clipped <- refl.scaled[, subset.index['xMin']:subset.index['xMax'],
                            subset.index['yMin']:subset.index['yMax']]



# plot the subset image RGB composite
rgb_ras <- plot_hs_rgb(refl = refl.clipped, 
            wavelengths = wavelengths,
            rgb.bands = rgb.bands, 
            proj4 = crs.info$Proj4, 
            ext = clip.extent,
            plt=TRUE)

# add the tree locations onto the RGB image
plotRGB(rgb_ras,
        r = 1, g = 2, b = 3,
        stretch = "lin",
        axes = TRUE,
        main="RGB Image",
        cex.main=2,
        addfun = points(polygons.in$X, 
                        polygons.in$Y, 
                        pch = 1, 
                        cex = polygons.in$crownDm, 
                        col = "white"))


# clip raster for single point
# get the index of the pixel in the image 
polygons.in[c('R', 'G', 'B')] <- raster::extract(x = rgb_ras,
                                                 y = cbind(polygons.in$X, polygons.in$Y))
  

# combine the data into a data frame to plot the spectral reflectances 
rgb.spectra.plotting <- polygons.in %>% 
  dplyr::select(indvdID, R, G, B) %>% 
  as.data.frame()

refl.gather <- tidyr::gather(rgb.spectra.plotting, 
                      key = "band",
                      value = "reflectance",
                      -indvdID) 

# add species information for each tree 
refl.gather <- merge(x = refl.gather, 
                             y = polygons.in[ ,c("indvdID","scntfcN")],
                             by = "indvdID") %>% 
  as.data.frame() %>% 
  dplyr::select(indvdID, scntfcN, everything())

# add wavelength values 
wavelength.lut <- as.data.frame(cbind(band = c('R','G','B'), wavelengths = rgb.bands)) 
refl.gather$wavelength <- wavelength.lut$wavelengths[match(refl.gather$band, wavelength.lut$band)]

# look at the first rows of this data frame 
head(refl.gather)

# plot the reflectance spectra 
ggplot(data = refl.gather, aes(x = wavelength, y = reflectance, colour = scntfcN)) +
  geom_line(aes(group = indvdID)) 



# how to calculate the point index within the image???!?!?!?!?!?!?!??!!??!!??!?!?!?
# to select ALL bands, not just RGB from the raster brick

i <- 1
point.index = {}
h5rows = full.extent@ymax - full.extent@ymin
point.index$x <- round(polygons.in$X[i] - full.extent@xmin)
point.index$y <- round(h5rows - (polygons.in$Y[i] - full.extent@ymin))
point.index


refl.pixel <- raster::extract(x = raster(refl.clipped),
                              y = polygons.points[1,c("X","Y")])







# clip raster for single polygon 




# clip using polygons ------------------------------------------------------
# adapt the Python tutorial here? http://neondataskills.org/HDF5/neon-aop-hdf5-py 
# or this one http://neondataskills.org/HDF5/Plot-Hyperspectral-Pixel-Spectral-Profile-In-R/



# display spectra ---------------------------------------------------------


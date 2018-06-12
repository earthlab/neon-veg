
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

# load custom functions
source("get_hs_band.R")


# Define the file name to be opened
f <- '~/github/neon-veg/NIWO/hyperspectral_reflectance/NEON_D13_NIWO_DP3_450000_4432000_reflectance.h5'

# look at the HDF5 file structure 
h5.struct <- h5ls(f,
                  all=T)

# read coordinate reference data
crs.info <- h5read(f, 
                   "NIWO/Reflectance/Metadata/Coordinate_System")
# convert "UTM" to lowercase "utm" for proper usage later
crs.info$Proj4 <- chartr("UTM", "utm", crs.info$Proj4)

# get attributes for the Reflectance dataset
refl.info <- h5readAttributes(f,
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

# wavelengths 
wavelengths <- h5read(f,
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
tile.extent <- extent(x.min, x.max, y.min, y.max)

# read reflectance data
# 1 band
band = 1 
refl <- h5read(f,"/NIWO/Reflectance/Reflectance_Data",
                 index = list(band, 1:n.cols, 1:n.rows))
# all bands
refl <- h5read(f,"/NIWO/Reflectance/Reflectance_Data",
                index = list(1:n.bands, 1:n.cols, 1:n.rows)) 

# view and apply scale factor to convert integer values to reflectance [0,1]
# and data ignore value
scale.factor <- refl.info$Scale_Factor
data.ignore <- refl.info$Data_Ignore_Value
refl[refl == data.ignore] <- NA 
refl <- refl / scale.factor


# plot single band based on wavelength  -------------------------------------------------------

# specify wavelength
wl <- 669
get_hs_band(refl, 
            wl,
            crs.info$Proj4,
            tile.extent)

source("get_hs_band.R")

# plot RGB composite ----------------------------------------------------
# specify indices of the chosen red, green, and blue bands to display 
bands <- c(669, 549, 474) 
plot_hs_rgb(refl, bands)


# create a function to process each band, create raster stack, plotRGB 



# clip using tree stem location points ------------------------------------




# clip using polygons ------------------------------------------------------
# adapt the Python tutorial here? http://neondataskills.org/HDF5/neon-aop-hdf5-py 
# or this one http://neondataskills.org/HDF5/Plot-Hyperspectral-Pixel-Spectral-Profile-In-R/



# display spectra ---------------------------------------------------------


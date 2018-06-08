
# This code is based off of the following NEON tutorials:
# Working with Hyperspectral Data in HDf5 Format in R
# http://neondataskills.org/HDF5/Imaging-Spectroscopy-HDF5-In-R/ 

# use the code below to install the rhdf5 library if it's not already installed.
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

# r Load `raster` and `rhdf5` packages and read NIS data into R
library(raster)
library(rhdf5)
library(rgdal)

# Define the file name to be opened
f <- '~/github/neon-veg/NIWO/hyperspectral_reflectance/NEON_D13_NIWO_DP3_450000_4432000_reflectance.h5'
# look at the HDF5 file structure 
h5_struct <- h5ls(f,all=T)

# read coordinate reference data
crsInfo <- h5read(f, "NIWO/Reflectance/Metadata/Coordinate_System")

# get attributes for the Reflectance dataset
reflInfo <- h5readAttributes(f,"NIWO/Reflectance/Reflectance_Data")

# resolution --> the Dimension Labels are "Line, Sample, Wavelength"
# does this correspond to nRows first?????? in the NEON tutorial,
# it was stated that the h5 reader reads cols first instead of rows.
# this won't be a problem for square rasters, but are the ones on the 
# edges not square? 
nRows <- reflInfo$Dimensions[1]
nCols <- reflInfo$Dimensions[2]
nBands <- reflInfo$Dimensions[3]

# wavelengths 
wavelengths <- h5read(f,
                      "NIWO/Reflectance/Metadata/Spectral_Data/Wavelength")

# read reflectance data
# Dimensions are listed as nBands x nRows (?) x nCols (?) in reflInfo
band = 1 # read 1 band. For all bands, use 1:nBands 
refl1 <- h5read(f,"/NIWO/Reflectance/Reflectance_Data",
               index = list(1:nBands, 1:nCols, 1:nRows)) 




# clip using polygons? 

# display spectral profiles (reflectance vs. wavelength)











# read in the wavelength information from the HDF5 file
wavelengths<- h5read(f,"wavelength")

# note that we can grab the dimensions of the dataset from the attributes
# we can then use that information to slice out our band data
nRows <- reflInfo$row_col_band[1]
nCols <- reflInfo$row_col_band[2]
nBands <- reflInfo$row_col_band[3]

# The HDF5 read function reads data in the order: Cols, Rows and bands
# This is different from how R reads data (rows, columns, bands). We'll adjust for 
# this later. 

#Extract or "slice" data for band 34 from the HDF5 file into array
b34 <- h5read(f,"Reflectance",index=list(1:nCols,1:nRows,34))

#Convert from array to matrix
b34 <- b34[,,1]

# plot a log of the data
image(log(b34))

# Plot range of reflectance values as a histogram to view range
# and distribution of values.
hist(b34,breaks=40,col="darkmagenta")

# View values between 0 and 5000
hist(b34,breaks=40,col="darkmagenta",xlim = c(0, 5000))

# reflectance values range from 0 to 1. The data scale factor in the 
# metadata tells us to divide all reflectance values by 10,000. 
# A value of 5,000 equates to a reflectance value of 0.50. 
# Storing as integers compared to floating points (with decimal places)
# creates a smaller file. 

# there is a no data value in our raster - let's define it
myNoDataValue <- as.numeric(reflInfo$`data ignore value`)
myNoDataValue

# set all values greater than 15,000 to NA
b34[b34 == myNoDataValue] <- NA

# Transpose x and y values in order for our final image to plot properly
b34<-t(b34)
image(log(b34), main="Transposed image")


# Create a georeferenced raster 

# Populate the raster image extent value. 
# get the map info, split out elements
mapInfo<-h5read(f,"map info")

# Extract each element of the map info information 
# so we can extract the lower left hand corner coordinates.
mapInfo<-unlist(strsplit(mapInfo, ","))

# grab resolution of raster as an object
res <- spInfo$xscale

# grab the left side x coordinate (xMin)
xMin <- as.numeric(mapInfo[4]) 

# grab the top corner coordinate (yMax)
yMax <- as.numeric(mapInfo[5])

# Calculate the lower right hand corner to define the full extent of the 
# raster. To do this we need the number of columns and rows in the raster
# and the resolution of the raster.

# note that you need to multiple the columns and rows by the resolution of 
# the data to calculate the proper extent!
xMax <- (xMin + (ncol(b34))*res)
yMin <- (yMax - (nrow(b34))*res) 

# define the extent (left, right, top, bottom)
rasExt <- extent(xMin,xMax,yMin,yMax)

# b34r isn't defined in the tutorial....
b34r <- raster(b34)

#assign the spatial extent to the raster
extent(b34r) <- rasExt

#look at raster attributes
b34r

#Create the projection in as object
myCRS <- spInfo$projdef

#define final raster with projection info 
#note that capitalization will throw errors on a MAC.
#if UTM is all caps it might cause an error!
b34r <- raster(b34,
               crs=myCRS)
extent(b34r) <- rasExt

# display image
image(log(b34r), 
      xlab = "UTM Easting", 
      ylab = "UTM Northing",
      main = "Properly Positioned Raster")

#let's change the colors of our raster and adjust the zlims 
col=terrain.colors(25)
image(b34r,  
      xlab = "UTM Easting", 
      ylab = "UTM Northing",
      main= "Raster w Custom Colors",
      col=col, 
      zlim=c(0,3000))

# write the raster as a geotiff
writeRaster(b34r,
            file="band34.tif",
            format="GTiff",
            overwrite=TRUE)

#It's always good practice to close the H5 connection before moving on!
#close the H5 file
H5close()

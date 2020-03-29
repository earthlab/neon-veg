# Classify vegetation species using NEON data
Perform tree species classification using freely available data from the National Ecological Observatory Network (NEON) and Random Forest models.

## Getting started 

### Packages to install

- devtools
- neonUtilities
- geoNEON 
- dplyr
- stringr
- sp
- raster
- rgeos
- tools
- randomForest
- gridExtra


### Workflow

In RStudio, click the "neon_veg_classify_species.Rproj" file to open this R project.

Open "01-main.R". 

### 1. Main
Install packages and define parameters for the other scripts. 
Download NEON woody vegetation data from the API and load it directly into R. 
Within this main script, run each of the following scripts in sequence to reproduce the analysis:

### 2. Create tree features
A point and circular polygon feature is created based on each tree location and maximum crown diameter measurement. These features are saved as shapefiles within the "data/data_output" folder. 

### 3. Process tree features
The proposed "clipping workflow" is applied to the tree polygons generated using in-situ NEON woody vegetation structure data:

1. Identical multi-bole entries are identified and excluded from subsequent analysis to remove duplicated points or polygons present in the raw data set.

2. An area threshold is applied to remove any small trees area values less than the area of two hyperspectral pixels.

3. "Occluded” polygons, those which are shorter and completely within the boundaries of other polygons, are removed from subsequent analysis Since they likely cannot be observed from the airborne perspective. Remaining polygons were checked for overlap with neighboring polygons. For each overlapping pair of polygons, shorter polygons were clipped by taller ones. If the remaining clipped area was smaller than the aforementioned area threshold, it was deleted. 

At this point, this workflow has generated a collection of polygons that will theoretically intersect with independent pixels in the airborne remote sensing data.

![Image of Clipping Workflow](https://github.com/earthlab/neon-veg/blob/master/images/clipping%20workflow.png)

### 4. Download AOP imagery

Download selected NEON Airborne Observation Platform (AOP) remote sensing mosaic data product tiles for the site and year of interest (saved to the following directory in the project: "data/data_raw/"). Each data product is downloaded to a deeply nested subdirectory structure in a top folder named with the data product ID (i.e. "DP3.30006.001"). Each set of files are then moved into a folder with a short intuitive pathname (i.e. "data/data_raw/hyperspectral")). 

| Data Product Name                                    | Data Product ID |
| :---                                                   | ---             |
| Ecosystem Structure (Canopy Height Model)             | DP3.30015.001   |
| Slope and Aspect - LiDAR                              | DP3.30025.001   |
| High-resolution orthorectified camera imagery  | DP3.30010.001   |
| Vegetation Indices - spectrometer - mosaic            | DP3.30026.001   |
| Spectrometer orthorectified surface directional reflectance | DP3.30006.001| 

### 5. Prep AOP imagery

Read each of the remote sensing data products and stack them into a single data cube per tile. This includes an aggregation operation to bring the high-resolution camera imagery (10cm x 10cm grid) to the 1m x 1m grid as used by the other data products. 

### 6. Plot AOP imagery

For a specified AOP imagery tile, create plots to visualize the remote sensing data. 
These plots include a RGB composite using hyperspectral bands, canopy height model, high-resolution digital camera RGB image:

![Image of AOP images](https://github.com/earthlab/neon-veg/blob/master/images/aop%20images.png) 


### 7. Extract training features 

Extract features (remote sensing image data values) for each sample (pixel) within the specified shapefile (containing points or polygons that correspond to trees at the NEON site).


### 8. Classify species 

Train a Random Forest (RF) model to classify tree species using in-situ tree species data as labels and remote sensing data as descriptive features. 


### 9. Assess accuracy

Compare classification accuracies across a series of RF models and present them in a table. 
Create a confusion matrix to show the classification accuracy for each species, in addition to user's and producer's accuracy metrics. 

### 10. Create figures 

Create figures for the manuscript: location of the NEON NIWO site within the United States, along with ribbon plots to illustrate the spectral signatures extracted for each species. 

neon-veg
================

A collection of R codes to create crown polygons based on field-based [NEON Woody Vegetation Structure 
](http://data.neonscience.org/api/v0/documents/NEON_vegStructure_userGuide_vA) data.

Installation
================

### Dependencies

To run the workflow, you'll need the following R packages:

-   devtools
-	geoNEON
-	sp
-	swfscMisc
-	rgdal
-	dplyr
- 	rgeos

Individual tree stems are recorded using distance and azimuth from a subplot reference point. To derive stem location in UTM coordinates, pull geolocation data for the sampling plots using the [geoNEON](https://github.com/NEONScience/NEON-geolocation/tree/master/geoNEON)  package.

To install the geoNEON package: 

``` r
library(devtools)
install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
library(geoNEON)
```

Setup
================

1. Download the Woody Plant Vegetation Structure data. If you downloaded it as a .zip file, be sure to unzip the file before proceeding. 


2. Within the SETUP section of the `main.R` script: 

- 	Set your current working directory to neon-veg in your local environment so the local functions can be loaded. 

```{r}
setwd("path/neon-veg")
```

-	Define the directory where the Woody Plant Vegetation Structure data sits as the `main_path` variable. 

-	Specify the output directory as the `out_dir` variable.

3. Run the `main.R` script to generate the following files in the output directory: 

- 	vst_merged.csv 
	*data table containing the merged vegetation structure entries based on individual ID)*

-	list_tiles.txt (a list of the 1km x 1km tiles containing woody veg stems, corresponds to the Canopy Height Model and RGB mosaic NEON data products)
-	polygons.shp (polygons generated with accompanying attribute (.dbf), shape index (.shx), and projection (.prj) files) 
-	polygons_checked_overlap.shp (previous shapefile after checking polygons for overlap and removing/editing them as necessary)
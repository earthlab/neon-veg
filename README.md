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

1. Set your current working directory to neon-veg in your local environment so the local functions can be loaded: 

```{r}
setwd("path/neon-veg")
```


2. After downloading the Woody plant vegetation structure data from the NEON Data Portal as a .zip file, define the directory where it sits as the `main_path` variable within the SETUP section of the main.R script. 


3. Specify the output path and filename for the shapefile to be generated as the `shp_filename` variable within the SETUP section of the main.R script. 
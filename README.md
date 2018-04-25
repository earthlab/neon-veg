neon-veg
================

A collection of R code to create crown polygons based on field-based NEON Woody Vegetation Structure data.


Installation
================

### Dependencies

To run the workflow, you'll need the following R packages:

- devtools
-	geoNEON
-	sp
-	swfscMisc
-	rgdal
-	dplyr
-	rgeos

Individual tree stems are recorded using distance and azimuth from a subplot reference point. To derive stem location in UTM coordinates, pull geolocation data for the sampling plots using the [geoNEON](https://github.com/NEONScience/NEON-geolocation/tree/master/geoNEON)  package.

To install the geoNEON package: 

``` r
devtools::install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
```

Setup
================

1. Download the Woody Plant Vegetation Structure data. If you downloaded it as a .zip file, be sure to unzip the contents into the directory for this project before proceeding (e.g., `neon-veg/SJER/`). 

*For a site with multiple dates of field data collection, there should be a series of folders with the collection date included in their names. Inside each folder, there are a series of Excel data tables including **vst_apparentindividual**, **vst_mappingandtagging**, and **vst_plotperyear**, which contain crown structure measurements, stem locations + species, and projection information, respectively. For more information, consult the [NEON User Guide to Woody Plant Vegetation Structure](http://data.neonscience.org/api/v0/documents/NEON_vegStructure_userGuide_vA).* 

We will assume that the working directory is set to the location of the `neon-veg` project.


Output
================

Run the `main.R` script to generate the following files in the output directory: 

-	**vst_merged.csv** - data table containing the merged vegetation structure entries based on individual ID

-	**list_tiles.txt** - a list of the 1km x 1km tiles containing woody veg stems, corresponds to the Canopy Height Model and RGB mosaic NEON data products

-	**polygons.shp** - polygons generated with accompanying attribute (.dbf), shape index (.shx), and projection (.prj) files

-	**polygons_checked_overlap.shp** - previous shapefile after checking polygons for overlap and removing/editing them as necessary

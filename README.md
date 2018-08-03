neon-veg
================

A collection of R code to create crown polygons based on field-based NEON Woody Vegetation Structure data.
A circular polygon is created for each tree with a mapped stem location, height, and crown diameter. 


Installation
================

### Dependencies

To run the workflow, you'll need the following R packages:

-   devtools
-	geoNEON
-   neonUtilities
-	sp
-	swfscMisc
-	rgdal
-	dplyr
-	rgeos
-	stringr
-	ggplot2
-	tidyr
-	dplyr
-	purrr
-	broom

The Woody Vegetation Structure data product is delivered in a series of data tables. The "mapping and tagging" table contains tree stem locations. The "apparent individual" table contains structural measurements of individual plants. The "plot per year" data table contains useful information about that year's sampling. For more information about these data, consult the [NEON User Guide to Woody Plant Vegetation Structure](http://data.neonscience.org/api/v0/documents/NEON_vegStructure_userGuide_vA). For a NEON site with multiple dates of field data collection, a series of data tables is generated for each date of collection. To combine these data into a single set of data tables, NEON provides the *stackByTable* function within the [neonUtilities](https://github.com/NEONScience/NEON-utilities/tree/master/neonUtilities) package.

To install the neonUtilities package: 

``` r
devtools::install_github("NEONScience/NEON-utilities/neonUtilities", dependencies=TRUE)
```

Individual tree stems are recorded using distance and azimuth from a subplot reference point. To derive stem location in UTM coordinates, pull geolocation data for the sampling plots by querying the NEON API using the [geoNEON](https://github.com/NEONScience/NEON-geolocation/tree/master/geoNEON) package.

To install the geoNEON package: 

``` r
devtools::install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
```


Setup
================

1. Download the Woody Plant Vegetation Structure data for your site and time period of interest. When downloaded from the [NEON Data Portal](http://data.neonscience.org/home), this data product is named **"NEON_struct-woody-plant.zip"** by default. 
2. Navigate to the location of the `neon-veg` project on your computer. 
3. Place the .zip veg structure data file inside the `neon-veg/data/` directory. 
4. Open the `main.R` script and edit the `data_path` and `out_dir` variables to reflect your file organization: 
-   `data_path` is the path and filename of the veg structure input data, such as *"data/NEON_struct-woody-plant.zip"*
-   `out_dir` is where the output files will be written. 

We will assume that the working directory is set to the location of the `neon-veg` project.


Output
================

Run the `main.R` script to generate the following files in the output directory: 

-	**vst_merged.csv** - data table containing the merged vegetation structure and mapping entries based on individual ID

-	**list_tiles.txt** - a list of the 1km x 1km tiles containing woody veg stems, corresponds to the Canopy Height Model and RGB mosaic NEON data products

-	**tree_counts.txt** - number of trees left after each major step in the workflow 

Five shapefiles generated with accompanying attribute (.dbf), shape index (.shx), and projection (.prj) files:

-	**mapped_stems.shp** - points for every mapped stem in the data set

-   **polygons_all.shp** - polygons for all complete entries (trees with mapped stem location, species, height, and crown diameter)

-   **polygons_filtered.shp** - previous shapefile after applying an area threshold (currently crowns with an estimated area of at least 4 hyperspectral pixels)

-	**polygons_checked_overlap.shp** - previous shapefile after checking polygons for overlap and removing/editing them as necessary

-	**mapped_stems_final.shp** - tree stem points for the final polygons after checking overlap


# create a directory for figures 
dir.create("figures")

library(data.table)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(USAboundaries)
library(ggspatial)


# Study Area ---------------------------------------------------------------

# NEON plot boundary for the NIWO site
# downloaded from: https://www.neonscience.org/data/about-data/spatial-data-maps
terrestrial_bounds <- sf::st_read("data/data_raw/neon_spatial_data_maps/fieldSamplingBoundaries/terrestrialSamplingBoundaries.shp")
niwo_bounds <- terrestrial_bounds[terrestrial_bounds$siteName == "Niwot Ridge Mountain Research Station",]

# NEON domain boundaries 
neon_domains <- sf::st_read("data/data_raw/neon_spatial_data_maps/NEONDomains_0/NEON_Domains.shp")

# create a SF point for the location of the NIWO site 
# latitude: 40.05425, longitude: -105.58237
niwo_point <- data.table(
  place=c("NIWO"),
  longitude=c(-105.58237),
  latitude=c(40.05425))

# convert to SF format for easy plotting with ggplot 
niwo_point_sf <- sf::st_as_sf(niwo_point, 
                              coords = c("longitude", "latitude"), 
                              crs = "+proj=longlat +datum=WGS84")


# Tutorial about maps in R: https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html

# get country data for the entire world 
world <- ne_countries(scale = "medium", returnclass = "sf")

# get states data (admin level 1 in USA)
states <- USAboundaries::us_states()

# isolate a single state from the bunch
colorado <- states[states$name == "Colorado",]

# make a map with a scale bar and north arrow using the ggspatial package 
ggplot() +
  geom_sf(data = world) +
  geom_sf(data = states) + 
  geom_sf(data = colorado, color = "black", fill = "darkgrey") + 
  geom_sf(data = niwo_point_sf, color = "black", size = 2) + 
  #geom_sf(data = niwo_bounds) + 
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-125, -70), ylim = c(25, 50)) + 
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() 

# write map to an image file 
ggsave(filename = file.path(dir_results, "study_area.png"), 
       device = "png", width = 8, height = 5)



# Ribbon plots ------------------------------------------------------------

# Ribbon plots 

# loop through each shapefile name, read the .csv containing spectral reflectance 
# for all trees within each data set, generate a ribbon plot and write to image file
# in the figures/ output directory 
for(i in 1:length(dirs_to_assess)){
  shapefile_filename <- file.path(dir_data_out, dirs_to_assess[i])
  print(paste("creating ribbon plots for:", dirs_to_assess[i]))
  
  createRibbonPlot(wavelengths, shapefile_filename, dir_data_out)
  
  createSeparateRibbonPlots(wavelengths, shapefile_filename, dir_data_out)
  
}


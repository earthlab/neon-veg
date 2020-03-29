# In this script, processing steps are applied to the tree features 
# generated using in-situ NEON woody vegetation structure data: 
# 
# 1. Identical multi-bole entries are identified and excluded 
# from subsequent analysis to remove duplicated points or polygons present 
# in the raw data set. 
# 
# 2. An area threshold is then applied to remove any small trees area values 
# less than the area of four hyperspectral pixels. This threshold 
# was selected with the coarser resolution of the hyperspectral and LiDAR data 
# products in mind. By preserving the trees with larger crowns, it is believed 
# that purer spectra will be extracted for them for the training data sets, 
# as opposed to extracting mixed pixels which signal from smaller plants and 
# background materials or neighboring vegetation. 
# 
# 3. “Engulfed” polygons, those which are shorter and completely within the 
# boundaries of other polygons, are also present in the initial crown polygons.
# Since they likely cannot be observed from the airborne perspective, 
# “engulfed” polygons were removed from subsequent analysis. Remaining 
# polygons were checked for overlap with neighboring polygons. For each 
# overlapping pair of polygons, shorter polygons were clipped by taller 
# ones. If the remaining clipped area was smaller than the aforementioned area 
# threshold, it was deleted. 
#
# At this point, this workflow has generated a collection of polygons that 
# will theoretically intersect with independent pixels in the airborne remote 
# sensing data.


# Remove multi-bole entries -----------------------------------------------

# Identify multi-bole entries: sometimes, there are multiple entries with 
# identical coordinates, height, and crown diameter (but their individualID's 
# are different, with "A", "B", etc. appended on the end of the last five 
# numbers.) In this step, the individualID values are assessed to find 
# multi-bole sets. If the coordinates, height, and crown diameters are 
# identical, then the entry with a letter appended on the end is deleted 
# from the analysis. This prevents duplicate polygons being used to extract 
# spectra. 

# Remove multibole entries from the merged woody veg data frame 
veg_multibole_removed <- remove_multibole(df = veg_merged)


# Count how many trees are left after removing multi-bole entries
tree_counts <- rbind(tree_counts
                    ,data.frame(count = c(nrow(veg_multibole_removed))
      ,description = c(" entries remain after removing multi-bole entries")))


# Convert the data frame to an SF object 
veg_multibole_removed_sf <- sf::st_as_sf(veg_multibole_removed
                                         ,coords = c("easting", "northing")
                                         ,crs = coord_ref)

# Create circular polygons with diamter = maximum crown diameter 
polys_multibole_rmvd_max_diam <- sf::st_buffer(veg_multibole_removed_sf
            ,dist = round((veg_multibole_removed_sf$maxCrownDiameter/2)
                                                          ,digits = 1))

# Create circular polygons with diamter = half of the maximum crown diameter 
polys_multibole_rmvd_half_diam <- sf::st_buffer(veg_multibole_removed_sf
            ,dist = round((veg_multibole_removed_sf$maxCrownDiameter/4)
                                                             ,digits = 1))


# Apply area threshold to remove small polygons ---------------------------

# Define the area threshold in units of [m^2]
# Area of RGB pixel, 10cm by 10cm
px_area_rgb <- 0.1 * 0.1 #[m^2]
# Gridded LiDAR products and HS pixels are 1m x 1m
px_area_hs <- 100 * px_area_rgb #[m^2]
# Multiply area of 1 HS pixel by the number of pixels; 
# as defined in the main script by the px_thresh variable
thresh <- px_area_hs * px_thresh #[m^2]

# Remove max-diam polygons with area < n hyperspectral pixels 
polygons_thresh_max_diam <- polys_multibole_rmvd_max_diam %>% 
  dplyr::mutate(area_m2 = sf::st_area(polys_multibole_rmvd_max_diam)) %>% 
  dplyr::filter(as.numeric(area_m2) > thresh)

# Remove half-diam polygons with area < n hyperspectral pixels 
polygons_thresh_half_diam <- polys_multibole_rmvd_half_diam %>% 
  dplyr::mutate(area_m2 = sf::st_area(polys_multibole_rmvd_half_diam)) %>% 
  dplyr::filter(as.numeric(area_m2) > thresh)


# Clip shorter polygons with taller polygons -----------------------------

polygons_clipped_max_diam <- clip_overlap(polygons_thresh_max_diam, thresh)
polygons_clipped_half_diam <- clip_overlap(polygons_thresh_half_diam, thresh)

# Check and fix/delete invalid geometries. 
# Remove invalid geometries if the reason
# for invalidity is "Too few points in geometry component[453729.741 4433259.265]"
# since there are too few points to create a valid polygon. 
# Use the reason = TRUE paramater in sf::st_isvalid to see the reason. 
polygons_clipped_max_diam_is_valid <- sf::st_is_valid(x = polygons_clipped_max_diam)
polygons_clipped_half_diam_is_valid <- sf::st_is_valid(x = polygons_clipped_half_diam)

# polygons_clipped_valid <- lwgeom::st_make_valid(polygons_clipped)
polygons_clipped_max_diam_valid <- polygons_clipped_max_diam %>% 
  dplyr::filter(polygons_clipped_max_diam_is_valid)
polygons_clipped_half_diam_valid <- polygons_clipped_half_diam %>% 
  dplyr::filter(polygons_clipped_half_diam_is_valid)

# VS-NOTE: need to test; run through the clipping workflow once more 
#polygons_clipped_max_diam2 <- clip_overlap(polygons_clipped_max_diam_valid, thresh)
#polygons_clipped_half_diam2 <- clip_overlap(polygons_clipped_half_diam_valid, thresh)


# Count how many trees are left after clipping areas of overlap and applying 
# area threshold to clipped polygons
tree_counts <- rbind(tree_counts
                     ,data.frame(count = c(nrow(polygons_clipped_max_diam_valid))
     ,description = c(" max_diam polygons remain after clipping overlap regions")))
tree_counts <- rbind(tree_counts
                     ,data.frame(count = c(nrow(polygons_clipped_half_diam_valid))
     ,description = c(" half_diam polygons remain after clipping overlap regions")))

# Write shapefile with clipped tree crown polygons
# maximum crown diameter
sf::st_write(obj = polygons_clipped_max_diam_valid
             ,dsn = file.path(dir_data_out, "veg_polys_max_diam_clipped_overlap.shp")
             ,delete_dsn = TRUE)
# half the maximum crown diameter 
sf::st_write(obj = polygons_clipped_half_diam_valid
             ,dsn = file.path(dir_data_out, "veg_polys_half_diam_clipped_overlap.shp")
             ,delete_dsn = TRUE)

# Write shapefile with a POINT for every mapped stem. 
# These points correspond to the half-diameter polygons. 
# Filter down the veg entries for those in the clipped half diam polygon set.
stems <- veg_multibole_removed %>% 
  filter(individualID %in% polygons_clipped_half_diam_valid$individualID)
# Create point geometries
points <- sf::st_as_sf(x = stems
                       ,coords = c("easting", "northing")
                       ,crs = coord_ref)
sf::st_write(obj = points
             ,dsn = file.path(dir_data_out, "veg_points_half_diam_clipped_overlap.shp")
              ,delete_dsn = TRUE)

# Write tree_counts information to a text file to assess the number of trees
# left after each processing step
write.csv(tree_counts, file = file.path(dir_data_out, "tree_counts.csv"))



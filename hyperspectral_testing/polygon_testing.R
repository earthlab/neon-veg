woody_merged <- merge(woody_utm, 
                      woody_ind_complete, 
                      by = "individualID") %>% 
  dplyr::filter(!is.na(height) & !is.na(maxCrownDiameter))


# create circular polygon for each stem based on max crown diameter 
woody_df_to_shp(df = woody_merged, 
                coord_ref = coord_ref,
                shrink = 1,
                num_sides = 24,
                shp_filename = paste(out_dir,
                                     "polygons_all",
                                     sep = ""))

# remove polygons with area < 4 hyperspectral pixels 
woody_thresh <- apply_area_threshold(woody_merged,
                                     nPix = 4)

# create circular polygon for each stem based on max crown diameter
woody_polygons <- woody_df_to_shp(df = woody_thresh, 
                                  coord_ref = coord_ref,
                                  shrink = 1,
                                  num_sides = 24,
                                  shp_filename = paste(out_dir,
                                                       "polygons_filtered",
                                                       sep = ""))

# multiply area of 1 pixel by the numPix input parameter
# rgb pixels are 25cm, 16 of them in a single HS pixel
thresh <- 25 * 25 * 16 * nPix # [cm^2]

# order polygons from tallest to shortest
df <- woody_polygons
polys_ordered <- df[order(df$height, 
                          decreasing = TRUE),]


# delete engulfed ---------------------------------------------------------
# find engulfed polygons (compare each pair of polygons. if area of overlap 
# is eqiuvalent to the area of either polygon, remove it from the list.)
polys_filtered <- polys_ordered

# create an empty list of polygon pairs that have been compared 
compared_pairs <- list();
c <- 1 # index for appending to list

for (individualID in polys_ordered@data$individualID){
  
  print(individualID)
  
  # if current polygon is not in filtered set, skip to the next current polygon
  if(sum(polys_filtered$individualID==individualID) == 0){
    next
  }
  
  # extract current vs. all other polygon from data frame
  current_poly <- get_poly(polys_filtered, index_type = 'id', number = individualID)
  other_polys <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
  
  # check for overlap between current polygon and all polygons
  overlap <- raster::intersect(current_poly, other_polys)
  n_overlap <- length(overlap)
  
  print(paste('n_overlap: ', as.character(n_overlap)))
  
  # if there is at least one polygon that overlaps with the current polygon,
  # check to see if it is engulfed. 
  if(n_overlap>0){
    
    # loop through the engulfed polygons. since the polygons have been reordered
    # from tallest to shortest, we can assume that any intersecting polygons 
    # are shorter (or of equal height). 
    
    for (o in 1:n_overlap){
      
      print(o)
      
      # check to see if the test polygon is still in the filtered polygon set... 
      test_poly <- get_poly(polys_filtered, 
                            index_type = 'id', 
                            number = overlap@data$individualID.2[[o]])
      # if test polygon is not in filtered set, skip to the next overlapping polygon
      if(sum(polys_filtered$individualID==test_poly$individualID) == 0){
        next
      }
      
      # get area of current and test polygons that overlap
      current_poly.area <- round(current_poly@polygons[[1]]@area, 2)
      test_poly.area <- round(test_poly@polygons[[1]]@area, 2) 
      
      print(paste("current_poly.area =", as.character(current_poly.area)))
      print(paste("test_poly.area =", as.character(test_poly.area)))
      
      # get area of the overlap between current polygon and test polygon
      overlap.area <- round(raster::intersect(current_poly, test_poly)@polygons[[1]]@area, 2)
      print(paste("overlap.area =", as.character(overlap.area)))
      
      clip <- rgeos::gIntersection(current_poly, test_poly, byid = TRUE, drop_lower_td = TRUE)
      
      # if total overlap
      if(current_poly.area == overlap.area | 
         test_poly.area == overlap.area){
        
        # then remove test polygon from filtered list
        polys_filtered <- polys_filtered[polys_filtered$individualID!=test_poly$individualID,]
        
        print("totally engulfed. removing test poly from polys_filtered")
        
      }
      
      }
    
  }
  
}

# write the filtered polygons to a shapefile. engulfed polygons should be gone. 
# write polygon(s) to shapefile  
suppressWarnings(
  writeOGR(polys_filtered, 
           getwd(),
           "polygons_engulfed_removed", 
           driver="ESRI Shapefile", 
           overwrite_layer = TRUE))


# clip shorter ------------------------------------------------------------
# for the remaining polygons, find the poly IDs that overlap
# clip shorter polygons; the taller polygons prevail
# add the pairs of polygon IDs to a list so they aren't compared again
# check area of clipped polygon(s) based on area threshold 


df_to_shp_points = function(df, coord_ref, shp_filename){
  # Creates a point shapefile for veg entries (rows) 
  # within the specified data frame of NEON vegetation structures. 
  #
  # Args: 
  #   df: data frame containing NEON veg structure entries that each
  #         include the columns: easting, northing, scientificname,
  #         taxonid, individualid
  #   coord_ref: Class "CRS" of Coordinate Reference System Arguments.
  #         Can be obtained from existing R layer with existing
  #         "coord .ref" field by using layer@crs
  #   shp_filename: string; output filename to write shapefile 
  # 
  # Returns:
  #   spdfs:
  #         spatialPointsDataFrame containing mapped stems as points
  
  message("\nCreating points based on mapped stem locations...")
  
  # select columns of interest
  stem_locations <- df %>%
    dplyr::select(easting, northing, individualID, scientificName, taxonID)
  
  # assign UTM coordinates to create SpatialPointsDataFrame
  coordinates(stem_locations) <- ~easting+northing
  proj4string(stem_locations) <- coord_ref
  
  # write polygon(s) to shapefile  
  suppressWarnings(
    writeOGR(stem_locations, 
             getwd(),
             shp_filename, 
             driver="ESRI Shapefile", 
             overwrite_layer = TRUE))
  
  return(stem_locations)
  
}

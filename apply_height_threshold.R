apply_height_threshold <- function(df, ht){ 
  # Applies a height threshold to remove polygons
  # with a height smaller than ht (meters)
  #
  # Args
  #   df
  #     data frame containing woody veg entries, including
  #     the height measurement
  #
  #   ht
  #     height (m) describing the height required to 
  #     keep polygons. 
  
  print("Removing polygons below height threshold...")
  
  # filter crowns with area < thresh
  df <- df %>%
          filter(height > ht) 
  
  return(df)
  
  }
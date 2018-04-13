polygon_overlap <- function(df, nPix){
  
  # multiply area of 1 pixel by the numPix input parameter
  # rgb pixels are 25cm, 16 of them in a single HS pixel
  thresh <- 25 * 25 * 16 * nPix # [cm^2]
  
  # order polygons from smallest to largest area
  polys_ordered <- df[order(df$maxCrownDiameter),]
  
  # create a polygon data frame to update with filtered/deleted/merged entries
  polys_filtered <- polys_ordered
  
  for (individualID in polys_ordered@data$individualID){
    
    # if this polygon was removed from the polys_filtered
    # data frame in a previous iteration,
    # skip it 
    if(sum(polys_filtered$individualID==individualID) == 0){
      break
    }
    
    # extract current polygon from data frame
    current_poly <- get_poly(polys_filtered, index_type = 'id', number = individualID)
    other_polys <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
    
    # check for overlap between current polygon and all polygons
    overlap <- intersect(current_poly, other_polys)
    n_overlap <- length(overlap)
    
    if(n_overlap>0){ 
      for (o in 1:n_overlap){
        
        # get area and height of current and test polygons that overlap
        current_poly.area <- round(current_poly@polygons[[1]]@area, 2)
        current_poly.height <- current_poly@data$height
        test_poly <- get_poly(polys_filtered, 
                              index_type = 'id', 
                              number = overlap@data$individualID.2[[o]])
        test_poly.area <- test_poly@polygons[[1]]@area
        test_poly.height <- test_poly@data$height
        
        # get area of the overlap between current polygon and test polygon
        overlap.area <- round(intersect(current_poly, test_poly)@polygons[[1]]@area, 2)
        
        # if total overlap (compare area of overlap to current polygon area),
        # remove current polygon from the collection of polygons
        if(current_poly.area == overlap.area | 
           test_poly.area == overlap.area | 
           ((current_poly.area - overlap.area) < 1)){
          polys_filtered <- other_polys
          break
        } 
        
        # partial overlap; different species
        else{ 
          # find overlapping area
          clip <- gIntersection(current_poly, test_poly, byid = TRUE, drop_lower_td = TRUE)
          
          # if the test polygon is taller, subtract the overlap area from the 
          # current polygon since it's shorter
          if(test_poly.height > current_poly.height){
            clipped <- (current_poly - clip)
            
            if(length(clipped) == 0){
              break
            }
            
            # if clipped area is large enough, save it in place of the current polygon 
            # and update the polys_filtered data frame 
            if((clipped@polygons[[1]]@area * 10000) > area_thresh){
              polys_filtered <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
              polys_filtered <- rbind(clipped, polys_filtered)
              # if the clipped area is too small, don't keep the current polygon in the 
              # polys_filtered data frame because it doesn't cover enough pixels! 
            } else{
              polys_filtered <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
            }
            
            # if the current polygon is taller, clip the test polygon
          } else if(test_poly.height > current_poly.height){
            clipped <- (test_poly - clip)
            
            # if clipped area is large enough, save it in place of the current polygon 
            # and update the polys_filtered data frame 
            # get individualID of test polygon
            j <- which(polys_filtered@data$individualID == test_poly@data$individualID)
            if((clipped@polygons[[1]]@area * 10000) > area_thresh){
              # replace with clipped polygon
              polys_filtered <- rbind(polys_filtered[1:j-1,],
                                      clipped, 
                                      polys_filtered[j+1:nrow(polys_filtered),])
              
              # if the clipped area is too small, don't keep the current polygon in the 
              # polys_filtered data frame because it doesn't cover enough pixels! 
            } else{
              polys_filtered <- polys_filtered[polys_filtered$individualID!=clipped@data$individualID,]
            }
          }
        }
      }
    }
  }
  writeOGR(polys_filtered, getwd(),
           "test_polygon_overlap", 
           driver="ESRI Shapefile", overwrite_layer = TRUE)

  
}
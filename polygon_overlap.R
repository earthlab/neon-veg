polygon_overlap <- function(df, nPix, shp_filename){
  
  message("\nChecking for overlap between polygons...")
  
  # multiply area of 1 pixel by the numPix input parameter
  # rgb pixels are 25cm, 16 of them in a single HS pixel
  thresh <- 25 * 25 * 16 * nPix # [cm^2]
  
  # order polygons from largest to smallest crown diameter
  polys_ordered <- df[order(df$crownDiam, 
                            decreasing = TRUE),]
  
  # create a polygon data frame to update with filtered/deleted/merged entries
  polys_filtered <- polys_ordered
  
  # create an empty list of polygon pairs that have been compared 
  compared_pairs <- list();
  c <- 1 # index for appending to list
  
  for (individualID in polys_ordered@data$individualID){
    
    # if this polygon was removed from the polys_filtered
    # data frame in a previous iteration, skip it 
    if(sum(polys_filtered$individualID==individualID) == 0){
      next
    }
    
    # extract current vs. all other polygon from data frame
    current_poly <- get_poly(polys_filtered, index_type = 'id', number = individualID)
    other_polys <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
    
    # check for overlap between current polygon and all polygons
    overlap <- raster::intersect(current_poly, other_polys)
    n_overlap <- length(overlap)
    
    if(n_overlap>0){ 
      for (o in 1:n_overlap){
        
        # if current polygon ID is not in filtered set
        if(sum(polys_filtered$individualID==individualID) == 0){
          break
        }
        
        # get area and height of current and test polygons that overlap
        current_poly.area <- round(current_poly@polygons[[1]]@area, 2)
        current_poly.height <- current_poly@data$height
        test_poly <- get_poly(polys_filtered, 
                              index_type = 'id', 
                              number = overlap@data$individualID.2[[o]])
        test_poly.area <- test_poly@polygons[[1]]@area
        test_poly.height <- test_poly@data$height
        
        # combine the ID's of the current and test polygons
        # to keep track of which pairs have been compared 
        id.pair <- paste(current_poly@data$individualID,
                         test_poly@data$individualID,
                         sep = " ")
        
        # if polygon pair was already compared, skip to the next pair
        if (id.pair %in% compared_pairs) {
          
          next
          
        } else { 
          
          # add to the list of polygon pairs that have been compared 
          compared_pairs[[c]] <- id.pair
          
          c <- c + 1 # increment index 
          
          # add opposite combination of polygons
          compared_pairs[[c]] <- paste(test_poly@data$individualID,
                                       current_poly@data$individualID,
                                       sep = " ")
          # increment index 
          c <- c + 1
        
          }
        
        # if test polygon is not in filtered set, skip to the next overlapping polygon
        if(sum(polys_filtered$individualID==test_poly$individualID) == 0){
          next
        }
        
        # get area of the overlap between current polygon and test polygon
        overlap.area <- round(raster::intersect(current_poly, test_poly)@polygons[[1]]@area, 2)
        clip <- gIntersection(current_poly, test_poly, byid = TRUE, drop_lower_td = TRUE)
        
        # if total overlap
        if(current_poly.area == overlap.area | 
                 test_poly.area == overlap.area | 
                 (current_poly.area - overlap.area) < 1){
                   if(current_poly.height > test_poly.height){ 
                      # then remove test polygon from filtered list
                      polys_filtered <- polys_filtered[polys_filtered$individualID!=test_poly$individualID,]
                   } else{ 
                     # remove current polygon
                     polys_filtered <- other_polys
                     break
                   }
          
                        
                 } else { # partial overlap
                   
                   # test polygon is taller
                   if(test_poly.height > current_poly.height){
                     clipped <- raster::erase(current_poly,
                                              raster::crop(current_poly, test_poly))

                   } else { 
                     # current polygon is taller 
                     clipped <- raster::erase(test_poly,
                                              raster::crop(test_poly, current_poly))
                   }
                   
                   # if there is no clipped area, skip to the next overlap polygon
                   if(length(clipped) == 0){
                     
                     next
                     
                   } else{ 
                     
                     # if clipped area exceeds area threshold
                     if(clipped@polygons[[1]]@area * 10000 > thresh){
                       
                       # replace current ploygon with clipped polygon
                       j <- which(polys_filtered@data$individualID == clipped$individualID)
                       
                       if(j==1){
                         polys_filtered <-  rbind(clipped, 
                                                   polys_filtered[(j+1):nrow(polys_filtered),])
                       
                         } else if(j==nrow(polys_filtered)) { # j is the last index

                           polys_filtered <- rbind(polys_filtered[1:(j-1),],
                                                   clipped)
                         
                       } else{
                         polys_filtered <- rbind(polys_filtered[1:(j-1),],
                                                 clipped, 
                                                 polys_filtered[(j+1):nrow(polys_filtered),])
                       }
                       
                       
                     } else{
                       # else, remove current polygon since clipped area is too small
                       polys_filtered <- polys_filtered[polys_filtered$individualID!=clipped$individualID,]
                     }
                     }
                 }
        }
    }
  }
  
print(polys_filtered)

# write final polygons to file after checking for overlap
writeOGR(polys_filtered, getwd(),
             paste(shp_filename),
             driver="ESRI Shapefile", overwrite_layer = TRUE)

return(polys_filtered)
  
}
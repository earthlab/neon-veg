get_poly = function(spdf, index_type, number){
  # this fuction extracts a single polygon inside a SpatialPolygonsDataFrame object
  # based on its individual ID OR index in the data frame for testing
  # index_type should be either "id" or "index"
  # if index_type == id, then it searches for the matching entry with that individual plant ID (number)
  # if index_type == index, then it uses the number as an index into the data frame
  
  if(index_type=='id'){ # use number as individual ID 
    i <- which(spdf$individualID == number)
  } else { # use number as index
    i <- number
  }
  
  print(i)
  
  coords = spdf@polygons[[i]]@Polygons[[1]]@coords
  extra_data = as.data.frame(spdf@data[spdf@data$ind_id == spdf$ind_id[i],], row.names = spdf$ind_id[i])
  
  # create SpatialPolygons
  P1 = Polygon(coords)
  Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = spdf$ind_id[i])), 
                        proj4string=spdf@CRS)
  
  # create SpatialPolygonsDataFrame
  Ps1 = SpatialPolygonsDataFrame(Ps1, 
                                 data = extra_data, match.ID = TRUE)
  return(Ps1)

  }


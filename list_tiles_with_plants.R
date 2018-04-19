list_tiles_with_plants <- function(woody,out_dir){
  # generate a list of 1km x 1km tiles containing field data
  # write list to a text file in the output directory provided
  
  print("Generating list of tiles containing stems...")
  
  # get easting, northing coordinates of all plants
  e <- NA
  n <- NA
  for (i in 1:dim(woody)[1]){
    easting <- floor(as.numeric(woody$easting[i])/1000)*1000
    northing <- floor(as.numeric(woody$northing[i])/1000)*1000
    e[i] <- easting
    n[i] <- northing
  }
  
  # find unique rows repesenting tiles
  easting_northings <- data.frame(as.character(e),as.character(n))
  colnames(easting_northings) <- c('e','n')
  tiles <- unique(easting_northings[,c('e','n')])
  
  # order by ascending tile coordinates 
  tiles <- tiles %>%
              arrange(e)
  
  # write to text file 
  tile_names <- paste(tiles$e, tiles$n, sep="_")
  tiles_file <- file(paste(out_dir,"list_tiles.txt", sep=""))
  writeLines(tile_names, tiles_file)
  close(tiles_file)
  
  return(tiles)
}
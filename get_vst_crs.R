get_vst_crs <- function(woody_path){
  
  # define path to vst_plotperyear table,
  # which contains geodetic datum and UTM zone 
  vst_path <- paste(woody_path, 
                      list.files(path = woody_path, 
                                 pattern = "plotperyear"), 
                      sep="/")
  
  vst_data <- read.csv(vst_path)
  
  datum <- as.character(vst_data$geodeticDatum[1])
  zone <- gsub("[^0-9\\.]", "", vst_data$utmZone[1])
  coord_ref <- CRS(paste("+proj=utm +zone=",zone, 
                         " +datum=",datum," +units=m",sep=""))
  return(coord_ref)
}
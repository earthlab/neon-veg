get_hs_band <- function(refl, wl, proj4, ext){
  
  # refl - array with reflectance data for all wavelengths
  # wl - wavelength(integer [nm]) to select
  # proj4 - projection string describing CRS for spatial object
  # ext - extent of the image to assign georeferences coordinates
  
  # get index of closest element in the wavelengths list
  wl.idx <- which.min(abs(wavelengths - wl))
  
  # extract reflectance data for the wavelength. convert array to matrix 
  refl.plot <- refl[wl.idx,,]
  
  # create raster and assign CRS
  refl.ras <- raster(refl.plot,
                    crs = proj4)
  
  # assign UTM extent to raster
  extent(refl.ras) <- ext
  
  print(refl.ras)
  
  # transpose x and y values for proper orientation in plot 
  refl.t <- t(refl.ras)
  
  # plot 
  image(log(refl.t), 
        main= paste("Band",
                    as.character(wl),
                    "nm",
                    sep= " "))
  
  return(refl.t)
  
  
}
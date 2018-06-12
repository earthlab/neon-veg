get_hs_band <- function(refl, wl, proj4, ext, plt){
  # extract a HS band with the specified wavelength (wl)
  # convert band to a raster, assign CRS and geographic extent, 
  # transpose for proper plotting orientation
  #
  # Args: 
  #   refl - array with reflectance data for all wavelengths
  #   wl - wavelength(integer [nm]) to select
  #   proj4 - projection string describing CRS for spatial object
  #   ext - extent of the image to assign georeferences coordinates
  #   plt - boolean variable to plot (TRUE) or not plot (FALSE) the band
  #
  # Returns:
  #   refl.t 
  #     transposed specified band of reflectance data with CRS and extent;
  #     ready to plot. 
  
  # get index of closest element in the wavelengths list
  wl.idx <- which.min(abs(wavelengths - wl))
  
  # extract reflectance data for the wavelength. convert array to matrix 
  refl.plot <- refl[wl.idx,,]
  
  # create raster and assign CRS
  refl.ras <- raster(refl.plot,
                    crs = proj4)
  
  # assign UTM extent to raster
  extent(refl.ras) <- ext
  
  # transpose x and y values for proper orientation in plot 
  refl.t <- t(refl.ras)
  
  # plot 
  if (plt == TRUE){
    image(log(refl.t), 
          main= paste("Band",
                      as.character(wl),
                      "nm",
                      sep= " "))
  }
  
  return(refl.t)
  
  
}
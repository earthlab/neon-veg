plot_hs_rgb <- function(refl, wavelengths, rgb.bands, proj4, ext, plt=FALSE){
  
  # extract the red, green, and blue bands from a hyperspectral 
  # data frame and plot the RGB image
  # 
  # Args: 
  #   refl
  #     array containing hyperspectral reflectance data 
  # 
  #   rgb.bands
  #     list of wavelengths (nm) for R, G, and B bands.
  #     for example: rgb.bands = c(620, 555, 450)
  
  # RED 
  r <- get_hs_band(refl = refl, 
                   wavelengths = wavelengths,
                   wl = rgb.bands[1],
                   proj4 = proj4,
                   ext = ext,
                   plt = FALSE)
  
  wl <- 549
  g <- get_hs_band(refl = refl, 
                   wavelengths = wavelengths,
                   wl = rgb.bands[2],
                   proj4,
                   ext,
                   plt = FALSE)
  
  wl <- 474
  b <- get_hs_band(refl = refl, 
                   wavelengths = wavelengths,
                   wl = rgb.bands[3],
                   proj4,
                   ext,
                   plt = FALSE)
  
  rgb_stack <- stack(r,g,b)
  rgb_brick <- brick(rgb_stack)
  
  if (plt == TRUE){
  plotRGB(rgb_brick,
          r = 1, g = 2, b = 3,
          stretch = "lin",
          axes = TRUE,
          main="RGB Image",
          cex.main=2)
  }
  
  return(rgb_brick)
  
  
}
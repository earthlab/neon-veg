# Plots AOP imagery for a specified tile.
# Optionally save the RGB composite images to file
# for the resambled high resolution digital camera imagery
# and the hyperspectral imagery. 

# specify which tile (.rds stacked AOP data file) to plot.
# 452000_4432000 has a visible road (defining feature) at the NIWO site
stacked_aop_data_filename = file.path(dir_data_out,
                                      "stacked_aop_data",
                                      paste0("stacked_aop_data_",
                                             east_north_string, ".rds"))

plot_aop_imagery <- function(filename, wavelengths, save_plots = FALSE){
  
  
  # read in the stacked AOP data for the specified tile
  stacked_aop_data <- readRDS(file = filename)
  
  # plot the high-res RGB image with mean values in each 1m by 1m cell
  rgb_stack <- raster::stack(stacked_aop_data$rgb_meanR, 
                             stacked_aop_data$rgb_meanG,
                             stacked_aop_data$rgb_meanB)
  rgb_brick <- raster::brick(rgb_stack)
  raster::plotRGB(rgb_brick,
                  r = 1, g = 2, b = 3,
                  stretch = "lin",
                  axes = TRUE,
                  main="Digital Camera RGB resampled to 1m grid",
                  cex.main=1)
  
  
  # plot RGB composite using selected bands from the hyperspectral imagery.
  # determine the index of each band closest to specific R,G,B wavelengths. 
  r_idx <- as.integer(which.min(abs(wavelengths - 620))) # Red 620 nm 
  g_idx <- as.integer(which.min(abs(wavelengths - 555))) # Green 555 nm 
  b_idx <- as.integer(which.min(abs(wavelengths - 450))) # Blue 450 nm 
  hs_rgb_stack <- raster::stack(stacked_aop_data@layers[[r_idx]],
                                stacked_aop_data@layers[[g_idx]],
                                stacked_aop_data@layers[[b_idx]])
  hs_rgb_brick <- raster::brick(hs_rgb_stack)
  
  raster::plotRGB(hs_rgb_brick,
                  r = 1, g = 2, b = 3,
                  stretch = "lin",
                  axes = TRUE,
                  main="Hyperspectral-derived \nRGB Composite using bands 620nm, 555nm, 450nm",
                  #xlab="Easting (m)",
                  #ylab="Northing (m)", 
                  cex.main=1)
  
  if(save_plots == TRUE){
    
    # write the resampled RGB to an image file for visualizing in QGIS
    writeRaster(rgb_brick,
                file.path(dir_data_out,"rgb_resampled_452000_4432000.tif"), 
                format="GTiff",
                overwrite=TRUE)
    
    # write the RGB composite to an image file for visualizing in QGIS
    writeRaster(hs_rgb_brick,
                file.path(dir_data_out,"hs_rgb_composite_452000_4432000.tif"), 
                format="GTiff",
                overwrite=TRUE)
    
  }
  
  
  
  # CHM 
  raster::plot(stacked_aop_data$chm,
               col=grey(1:100/100),
               axes = TRUE,
               main="Lidar-Derived Canopy Height Model",
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1, # title text size 
               legend.args=list(text='Height above ground [m]',side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar
  
  # slope
  raster::plot(stacked_aop_data$slope,
               col=grey(1:100/100),
               axes = TRUE,
               main="Lidar-Derived Slope",
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1,
               legend.args=list(text='Ratio of rise over run [degrees]',
                                side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar)
  
  # aspect 
  raster::plot(stacked_aop_data$aspect,
               col=grey(1:100/100),
               axes = TRUE,
               main="Lidar-Derived Aspect",
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1,
               legend.args=list(text='Direction of steepest slope \n [degrees from North]',
                                side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar)
  
  # NDVI vegetation index
  raster::plot(stacked_aop_data$NDVI,
               col=grey(1:100/100),
               axes = TRUE,
               main="Hyperspectral-Derived NDVI",
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1,
               legend.args=list(text='NDVI [0,1]',side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar))
  
}

plot_aop_imagery(stacked_aop_data_filename, wavelengths, save_plots = TRUE)

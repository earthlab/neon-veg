# This script trains a Random Forest (RF) model to classify tree species
# using the in-situ tree measurements for species labels 
# and remote sensing data as descriptive features. 


# get a description of the shapefile to use for naming outputs
shapefile_description <- tools::file_path_sans_ext(basename(shapefile_filename))
# define the csv file containing extracted features
extracted_features_filename <- file.path(dir_data_out,
                                         paste0(shapefile_description,
                                                "-extracted_features.csv"))


# description to name a folder for the classification outputs 
out_description <- paste0("rf_",shapefile_description) 
check_create_dir(file.path(dir_data_out,out_description))

# taxonIDs (species) to predict
taxon_list <- c("ABLAL","PICOL","PIEN","PIFL2")


# RF tuning parameter, number of trees to grow. default value 500
ntree <- 5000

# To reduce bias, this boolean variable indicates whether the same number 
# of samples are randomly selected per species class. Otherwise,
# all samples per class are used for training the classifier. 
randomMinSamples <- FALSE

# To remove the sample size bias, if TRUE this filters down each of the raw NEON 
# shapefile data sets to only contain the individualIDs present in the neon_veg set 
#neonvegIDsForBothShapefiles <- TRUE

# boolean variable. if TRUE, keep separate set for validation
#independentValidationSet <- TRUE 

# randomly select this amount of data for training, use the rest for validation
#percentTrain <- 0.8 

pcaInsteadOfWavelengths <- TRUE
nPCs <- 2 # number of PCAs to keep 

# keep most important variables and run RF again with reduced feature set 
#keepMostImpVar <- FALSE

# create boxplots and write to file 
#createBoxplots <- TRUE

# define the "bad bands" wavelength ranges in nanometers, where atmospheric 
# absorption creates unreliable reflectance values. 
bad_band_window_1 <- c(1340, 1445)
bad_band_window_2 <- c(1790, 1955)
# remove the bad bands from the list of wavelengths 
remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                               wavelengths < bad_band_window_1[2]) | 
                              (wavelengths > bad_band_window_2[1] & 
                                 wavelengths < bad_band_window_2[2])]

# create a LUT that matches actual wavelength values with the column names,
# X followed by the rounded wavelength values. Remove the rows that are 
# within thebad band ranges. 
wavelength_lut <- data.frame(wavelength = wavelengths,
                             xwavelength = paste0("X", 
                                                  as.character(round(wavelengths))),
                             stringsAsFactors = FALSE) %>% 
  filter(!wavelength %in% remove_bands) 

# features to use in the RF models.
# this list is used to filter the columns of the data frame,
# to remove the ones containing other metadata per tree from the model. 
featureNames <- c("taxonID", 
                  wavelength_lut$xwavelength,
                  "chm", 
                  "slope", 
                  "aspect",
                  "ARVI",
                  "EVI",
                  "NDLI",
                  "NDNI",
                  "NDVI",
                  "PRI",
                  "SAVI",
                  #"rgb_meanR",
                  #"rgb_meanG",
                  #"rgb_meanB",
                  "rgb_mean_sd_R",
                  "rgb_mean_sd_G",
                  "rgb_mean_sd_B",
                  "pixelNumber",
                  "eastingIDs",
                  "northingIDs")


# open a text file to record the output results
rf_output_file <- file(file.path(dir_data_out,out_description,
                                 "rf_model_summaries.txt"), "w")

# start the timer
start_time <- Sys.time()


print("Currently training the random forest with features extracted using:")
print(extracted_features_filename)

# read the values extracted from the data cube
df_orig <- read.csv(extracted_features_filename)

# filter down to the species of interest 
df_orig <- df_orig %>% 
  dplyr::filter(taxonID %in% taxon_list)

# # remove any pixels in the independent validation set
# if(independentValidationSet){
#   
#   # combine the pixelNumber, easting, and northing into a string
#   # for each row with underscore separation
#   dfIDs <- paste(df_orig[,"pixelNumber"],
#                  df_orig[,"eastingIDs"],
#                  df_orig[,"northingIDs"], 
#                  sep = "_") 
#   
#   # add this as a new column to the data frame 
#   df_orig <- df_orig %>% mutate(dfIDs = dfIDs)
#   
#   # remove any pixels that are in the independent validation set
#   df_orig <- df_orig[!(df_orig$dfIDs %in% valIDs),]
# 
# } 


# testing the influence of sampling bias - this code filters the "raw"
# NEON polygon sets down to contain the same individual IDs as the corresponding
# clipped datas set. This is meant to assess the potential influence of sampling
# bias (since there is nearly double # of samples for the raw data sets compared
# to those after the clipping steps), and see what effect the clipping steps 
# have on species classification using circular polygons
if(neonvegIDsForBothShapefiles){
  
  # check for the clipped half diam shapefile; keep track of the individual ID's 
  if(grepl("polygons_half_diam", extracted_features_filename)){
    print(paste0("Filtering individualIDs to match those in the clipped half diam polygons:", 
                 extracted_features_filename))
    
    # read in the clipped half diam polygons
    # VS-NOTE: make this a variable instead of hard-coding the file path and name 
    df_clipped <- read.csv(file.path(dir_data_out,
              "veg_polys_half_diam_clipped_overlap-extracted_features.csv"))
    
    # figure out which individual IDs are present in the clipped polygon data set
    clipped_IDs <- unique(droplevels(df_clipped$indvdID))
    filterIDs <- clipped_IDs
    
    # delete the clipped data frame from memory 
    rm(df_clipped)
    
    df_orig <- df_orig %>% dplyr::filter(indvdID %in% filterIDs) %>% droplevels()
  }
  
  
  # check for the clipped max diam shapefile; keep track of the individual ID's 
  if(grepl("polygons_max_diam", extracted_features_filename)){
    print(paste0("Filtering individualIDs to match those in the clipped max diam polygons:", 
                 extracted_features_filename))
    
    # read in the clipped max diam polygons
    # VS-NOTE: make this a variable instead of hard-coding the file path and name 
    df_clipped <- read.csv(file.path(dir_data_out,
                  "veg_polys_max_diam_clipped_overlap-extracted_features.csv"))
    
    # figure out which individual IDs are present in the clipped polygon data set
    clipped_IDs <- unique(droplevels(df_clipped$indvdID))
    filterIDs <- clipped_IDs
    
    # delete the clipped data frame from memory 
    rm(df_clipped)
    
    df_orig <- df_orig %>% dplyr::filter(indvdID %in% filterIDs) %>% droplevels()
  }
  
}



# Remove any spectra that have a height == 0
print(paste0(as.character(sum(df_orig$chm==0)), 
             " pixels have a height of 0 in the CHM"))
print("Removing these rows from the training set ... ")

# also reset the factor levels (in case there are dropped taxonID levels)
df <- df_orig %>% filter(chm>0) %>% droplevels()

# remove the bad bands from the list of wavelengths 
remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                               wavelengths < bad_band_window_1[2]) | 
                              (wavelengths > bad_band_window_2[1] & 
                                 wavelengths < bad_band_window_2[2])]

# create a LUT that matches actual wavelength values with the column names,
# X followed by the rounded wavelength values. Remove the rows that are 
# within thebad band ranges. 
wavelength_lut <- data.frame(wavelength = wavelengths,
                             xwavelength = paste0("X", 
                                                  as.character(round(wavelengths))),
                             stringsAsFactors = FALSE) %>% 
  filter(!wavelength %in% remove_bands) 

# filter the data to contain only the features of interest 
features <- df %>% 
  dplyr::select(featureNames)

# testing whether PCA yields better accuracy than individual wavelength reflectance data
if(pcaInsteadOfWavelengths == TRUE){
  
  # remove the individual spectral reflectance bands from the training data
  features <- features %>% dplyr::select(-c(wavelength_lut$xwavelength))
  
  print(colnames(features))
  
  # PCA: calculate Principal Components 
  hs <- df %>% dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()
  hs_pca <- stats::prcomp(hs, center = TRUE, scale. = TRUE)
  summary(hs_pca)
  features <- cbind(features, hs_pca$x[,1:nPCs]) # add first n PCs to features data frame
  # visualize where each sample falls on a plot with PC2 vs PC1 
  # ggbiplot::ggbiplot(hs_pca, 
  #                    choices = 1:2, # which PCs to plot
  #                    obs.scale = 1, var.scale = 1, # scale observations & variables 
  #                    var.axes=FALSE, # remove arrows 
  #                    groups = df$taxonID, # color the points by species
  #                    ellipse = TRUE, # draw ellipse around each group
  #                    circle = TRUE # draw circle around center of data set
  # )   + 
  #   ggtitle("PCA biplot, PC1 and PC2") + 
  #   scale_color_brewer(palette="Spectral") + 
  #   theme_bw()
  # # save to file 
  # ggplot2::ggsave(paste0(out_dir, outDescription, 
  #                        "pcaPlot_",shapefileLayerNames$description[i],
  #                        ".png"))
  
}

# VS-NOTE: need to scale feature values before training classifier
print("Features used in current RF model: ")
print(colnames(features))


# count the number of samples per species 
featureSummary <- features %>%
  dplyr::group_by(as.character(taxonID)) %>%
  dplyr::summarize(total = n()) 

print("number of samples per species class")
print(featureSummary)


if(independentValidationSet){
  
  # perform PCA
  # remove the individual band reflectances.
  # this only needs to be done once, so check if the validationSet
  # already has a column named "PC1". 
  if(pcaInsteadOfWavelengths == TRUE){
    message("in first loop")
    
    if(("PC1" %in% colnames(validationSet)) == FALSE){
      
      # filter the data to contain only the features of interest 
      validationSet <- validationSet %>% 
        dplyr::select(featureNames)
      
      print("performing PCA on validation set...")
      wl_removed <- validationSet %>% as.data.frame() %>% 
        dplyr::select(-c(wavelength_lut$xwavelength))
      
      # isolate the individual band reflectances
      val_hs <- validationSet %>% 
        dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()
      
      val_hs_pca <- stats::prcomp(val_hs, center = TRUE, scale. = TRUE)
      summary(val_hs_pca)
      validationSet <- cbind(wl_removed, val_hs_pca$x[,1:nPCs])
      
    }
    
  }
  
  
  
  # concatenate pixel #, easting, and northing 
  # for an individual entry per pixel in the input imagery,
  # i.e. "264231_452000_4431000"
  features$dfIDs <- paste(features[,"pixelNumber"],
                          features[,"eastingIDs"],
                          features[,"northingIDs"], 
                          sep = "_") 
  
  #print("The following pixelNumber_easting_northing values are ")
  #print("found in the current input set and the validation set: ")
  #print(intersect(features$dfIDs, valIDs))
  
  print("Originally this many rows in the features DF: ")
  print(nrow(features))
  
  # remove any spectra that are being kept separate for the 
  # independent validation set 
  features <- features[!(features$dfIDs %in% valIDs),]
  
  print("After removing pixels that are in the validation set, nrow = ")
  print(nrow(features))
  
  # remove the pixelNumber, easting, and northing columns since they
  # are not input features to the train the classifier 
  features <- features %>% 
    dplyr::select(-c(pixelNumber, eastingIDs, northingIDs, dfIDs))
  
}




if(randomMinSamples){
  # reduce number of samples per species to avoid classifier bias
  
  # count the minimum number of samples for a single class
  minSamples <- min(featureSummary$total)  
  print(paste0("Randomly selecting ",
               as.character(minSamples),
               " samples per species class to avoid classifier bias"))
  
  # isolate the samples per species
  taxon1 <- features[features$taxonID==taxon_list[1],]
  taxon2 <- features[features$taxonID==taxon_list[2],]
  taxon3 <- features[features$taxonID==taxon_list[3],]
  taxon4 <- features[features$taxonID==taxon_list[4],]
  
  # keep random minSamples of each species; merge
  taxon1 <- taxon1[sample(nrow(taxon1), minSamples), ]
  taxon2 <- taxon2[sample(nrow(taxon2), minSamples), ]
  taxon3 <- taxon3[sample(nrow(taxon3), minSamples), ]
  taxon4 <- taxon4[sample(nrow(taxon4), minSamples), ]
  
  features <- rbind(taxon1, taxon2, taxon3, taxon4)
  
} else{
  #print("Using all samples per class")
}



# TRAIN RF CLASSIFIER using training set  ---------------------------------
rf_startTime <- Sys.time()
set.seed(104)
rf_model <- randomForest::randomForest(as.factor(features$taxonID) ~ .,
                                       data=features, 
                                       importance=TRUE, 
                                       # ntree = # of trees to grow
                                       ntree=ntree) 
print("randomForest time elapsed for model training: ")
print(Sys.time()-rf_startTime)

print(rf_model)

# save RF model to file 
save(rf_model, file = file.path(dir_data_out, out_description,
                                paste0("rf_model_",shapefile_description,".RData")))

# write all relevant information to the textfile: 
# shapefile name
write(shapefile_description, rf_output_file, append=TRUE)
write("\n", rf_output_file, append=TRUE) #newline

# number of samples per class
featureSummary <- data.frame(featureSummary)
colnames(featureSummary) <- c("taxonID","numberOfSamples")
capture.output(featureSummary, file = rf_output_file, append=TRUE)

# y = predicted data; (horizontal axis)
# x = observed data (true class labels) (vertical axis)
accuracy <- rfUtilities::accuracy(x = rf_model$y,
                                  y = rf_model$predicted)

# record each accuracy metric in the table for a final comparison.
# round each value to the nearest decimal place 
OA_OOB <- round(accuracy$PCC, 1) # Overall Accuracy
K <- round(accuracy$kappa, 3) #Cohen's Kappa 


write("\nOverall Accuracy:", rf_output_file, append=TRUE) #newline
write(accuracy$PCC, rf_output_file, append=TRUE)

write("\nUser's Accuracy:", rf_output_file, append=TRUE) #newline
capture.output(accuracy$users.accuracy, file = rf_output_file, append=TRUE)

write("\nProducer's Accuracy:", rf_output_file, append=TRUE) #newline
capture.output(accuracy$producers.accuracy, file = rf_output_file, append=TRUE)

write("\nConfusion Matrix:", rf_output_file, append=TRUE) #newline
capture.output(accuracy$confusion, file = rf_output_file, append=TRUE)

write("\nCohen's Kappa:", rf_output_file, append=TRUE) #newline
capture.output(accuracy$kappa, file = rf_output_file, append=TRUE)


# INDEPENDENT VALIDATION  -------------------------------------------------

# predict species ID for validation set 
if(independentValidationSet){
  predValidation <- predict(rf_model, validationSet, type = "class")
  confusionTable <- table(predValidation, validationSet$taxonID)
  print(confusionTable)
  val_OA <- sum(predValidation == validationSet$taxonID) / 
    length(validationSet$taxonID)
  print(paste0("overall accuracy predicting validation set: ",
               as.character(val_OA)))
  # write the accuracy summary data frame to file 
  write.csv(confusionTable,
            file.path(dir_data_out, out_description, 
                      paste0("rfConfusionMatrix_independentValidationSet_",
                             shapefile_description,"_Accuracy_",
                             as.character(round(val_OA, 3)),".csv")))
  
}


# write all relevant information to the textfile: 
# features used to describe each sample (pixel)
write("\ndescriptive features used to train this model: ", rf_output_file, append=TRUE) #newline
write(colnames(features), rf_output_file, append=TRUE)

# RF model summary, OOB error rate 
capture.output(rf_model, file = rf_output_file, append=TRUE)

# close the text file
close(rf_output_file)

end_time <- Sys.time()
print("Elapsed time: ")
print(end_time-start_time)




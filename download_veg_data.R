# This script downloads NEON woody vegetation data from the API and loads 
# it directly into R. 


# Let's load the in-situ Woody Vegetation Structure data straight into R.
# Specify the NEON site and starting/ending date(s) for the data. 
# A message in the console will display the total file size to be downloaded.
# Proceed by typing "y" and pressing Enter. 
veg_raw <- neonUtilities::loadByProduct(dpID = "DP1.10098.001"   
                                        ,site = site_code              
                                        ,startdate = veg_start_date      
                                        #,enddate = "YYYY-MM"
                                        ,package = "basic"          
                                        ,check.size = T)

# ------------------------------------------------------------------------
# # Alternative option to loading data straight into R:
# # Download the files locally from NEON API using the zipsByProduct function. 
# # They will be stored within "data/data_raw/filesToStack#####" where "#####" 
# # is the middle portion of the data product ID string. 
# # Since the data product ID is "DP1.10098.001", files will be saved 
# # in "filesToStack10098". 
# neonUtilities::zipsByProduct(dpID = "DP1.10098.001"   
#                              ,site = "NIWO"              
#                              ,startdate = "2016-01"      
#                              #,enddate = "YYYY-MM"
#                              ,package = "basic"          
#                              ,savepath = "data/data_raw"
#                              ,check.size = T)
# 
# # The downloaded files can now be passed to stackByTable() to be stacked. 
# # They are saved within a folder called "stackedFiles" within the "filesToStack10098) directory. 
# neonUtilities::stackByTable(filepath = "data/data_raw/filesToStack10098" 
#                             ,folder = T)
# ------------------------------------------------------------------------
merge_vst_tables <- function(vst_mapping, vst_individual){
    # Merges species data with vegetation structure data
    # based on matching individual ID's. 
    # 
    # inputs
    #   vst_mapping: 
    #     data frame containing scientific names of mapped plants (vst_mappingandtagging)
    #   vst_individual:
    #     data frame containing vegetation structure (vst_apparentindividual)
    #
    # output
    #   merged:
    #     input data frame with additional columns of information
    
    
    # find matching individual ID's between the tables
    i <- match(as.character(vst_mapping$individualID),
               as.character(vst_individual$individualID))
    
    # check for duplicates and keep most recent measurement? 
    
    # add height and crown diameter columns to mapping table
    merged <- vst_mapping
    merged$height <- vst_individual$height[i]
    merged$maxCrownDiameter <- vst_individual$maxCrownDiameter[i]
    # replace the date column with apparentindividual 
    # since the location of the plant will not change but the height 
    # may between measurements
    merged$date <- vst_individual$date[i]
    
    # keep only entries that have height and crown diameter values
    merged <- merged[complete.cases(merged$height) & 
                     complete.cases(merged$maxCrownDiameter),]
    
    return(merged)
}
  

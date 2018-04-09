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
    
    # add species columns to structure table
    merged <- vst_mapping
    merged$height <- vst_individual$height[i]
    merged$maxCrownDiameter <- vst_individual$maxCrownDiameter[i]
    
    return(merged)
}
  

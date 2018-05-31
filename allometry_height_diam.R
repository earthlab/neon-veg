allometry_height_diam <- function(df){
  
  # calculates linear regression models for each taxonID 
  # within the input data frame. 
  # Dependent variable: crown diameter (m) 
  # Independent variable: height (m)
  
 # libraries: ggplot2, tidyr, purrr, dplyr, broom 
  

  # nest the data by taxon ID: 
  # create separate list of data for each species
  by_taxonID <- stems_final %>% 
    group_by(taxonID) %>% 
    nest()
  
  # function to calculate linear model of crownDiam 
  # as a function of height 
  taxon_model <- function(df){
    lm(crownDiam ~ height, data = df)
  }
  
  # for each species data set, calculate linear model
  models <- by_taxonID %>% 
    mutate(
      model = map(data, taxon_model)
    )
  
  # compute summary, r-squared, parameters, and observation statistics 
  # for each linear model using the "broom" library 
  models <- models %>% 
    mutate( 
      glance  = map(model, broom::glance),
      rsq     = glance %>% map_dbl("r.squared"),
      tidy    = map(model, broom::tidy),
      augment = map(model, broom::augment)
      )
  
  # plot data points and fitted models using ggplot 
  g <- models %>%
    unnest(data) %>%
    ggplot(aes(height, crownDiam)) +
    geom_point(shape=1) + 
    facet_wrap(~taxonID) +    
    geom_smooth(method=lm) +
    labs(x = "height (m)", y = "crown diameter (m)")
  print(g)
  
  return(models)
    
}

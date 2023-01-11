perc_agree_on_sign_weighted_AR4_method <- function(all_data, weights) {
  
  #Test if models agree on the sign of the change
  
  #"Robust" when magnitude of the multi-model ensemble mean exceeds the inter-model standard deviation
  
  #Initialise 
  agr <- NA
  
  #If pixel has all models available
  if (all (!is.na(all_data))) {
    
    #Set zero for non-robust)
    agr <- 0
    
    #Calculate weighted multi-model mean
    multimodel_mean <- weighted.mean(all_data, w=weights)
    
    
    #Calculate weighted inter-model standard deviation
    stdev <- wtd_stdev(all_data, weights=weights)
    
    
    if (abs(multimodel_mean) > stdev) agr <- 1
    
  }
  
  
  return(agr)
  
}



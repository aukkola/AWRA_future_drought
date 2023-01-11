perc_agree_with_obs_mean_weighted <- function(all_data, weights=NA, threshold=10) {
  
  #NB. Assumes first layer is mean_data, and the rest all_data !!!!!!!!!!
  mean_data <- all_data[1]
  all_data  <- all_data[2:length(all_data)]
  
  
  #If missing, return NA
  if (is.na(mean_data)) {
    return(NA)
  }
  
  
  #Calculate percentage error
  mean_error <- abs((all_data / mean_data) -1) * 100
  
  ind <- which(mean_error <= threshold)
  
  
  #Calculate agreement
  if (all(is.na(weights))) {
    
    agr <- length(ind)/length(all_data)
    
  } else {
    
    #Number of independent models
    no_independent <- sum(weights)
    
    #Calculate agreement from weights
    agr <- sum(weights[ind]) / no_independent
    
  }
  
  
  
  return(agr)
  
}
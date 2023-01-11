perc_agree_on_mean_weighted <- function(all_data, weights=NA, threshold=10) {
  
  
  
  #Calculate mean
  if (all(is.na(weights))) {
    mean_data <- mean(all_data, na.rm=TRUE)
    
  } else {
    mean_data <- weighted.mean(all_data, w=weights, na.rm=TRUE)
  }
  
  
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
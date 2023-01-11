perc_agree_on_sign_weighted <- function(all_data, weights) {
  
  #Test if models agree on the sign of the change
  
  
  #Get the sign of change
  sign_change <- sign(all_data)
  
  
  #Calculate which sign is produced by the majority of models
  
  half <- length(all_data) / 2
  
  #Majority positive
  if (length(which(sign_change > 0)) > half) {
    
    ind <- length(which(sign_change > 0))
    
  #Majority negative
  } else if (length(which(sign_change < 0)) > half) {
  
    ind <- length(which(sign_change < 0))
    
  #Equal positive and negative
  } else {
    
    #Set index as half of models so agreement becomes 50%
    ind <- length(all_data)/2
    
  }


  
  #Not weighted
  if (all(is.na(weights))) {
    
    agr <- ind / length(all_data)
    
  #Weighted
  } else {
    
    #Number of independent models
    no_independent <- sum(weights)
    
    #Calculate agreement from weights
    agr <- sum(weights[ind]) / no_independent
    
    
  }
  
  
  return(agr)
  
}



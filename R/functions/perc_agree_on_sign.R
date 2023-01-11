perc_agree_on_sign <- function(all_data) {
  
  #Calculate median
  median_data <- median(all_data, na.rm=TRUE)
  
  
  #If missing, return NA
  if (is.na(median_data)) {
    return(NA)
  }
  
  
  
  #If positive
  if (median_data > 0) {
    
    ind <- which(all_data > 0)
    
    return(length(ind)/length(all_data))
    
  #If negative
  } else if (median_data <= 0) {
    
    ind <- which(all_data < 0)
    
    return(length(ind)/length(all_data))
    
  }

}
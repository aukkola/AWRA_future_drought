#Calculate frequency
calc_freq <- function(x) {
  
  if (all(is.na(x))) {
    
    return(NA)
    
  } else {
    
    #Number of droughts events divided by total number of months
    freq <- length(which(!is.na(x))) / length(x)
    
    return(freq)
  }
  
}
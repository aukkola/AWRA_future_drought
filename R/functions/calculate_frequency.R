calculate_frequency <- function(data, window=120) {
  
  #Get data length
  len <- length(data)
  
  #Calculate number of time blocks
  start_ind <- seq(1, by=window, length.out=ceiling(len/window))
  
  end_ind <- seq(window, by=window, length.out=length(start_ind))
  
  #Adjust end index if longer than data
  if (max(end_ind) > len) end_ind[length(end_ind)] <- len

  
  #Acccount for possibly different length of last time block by
  #normalising with time block length
  scale <- end_ind - start_ind + 1
  
  
  #Calculate number of events per time block
  frequency <- mapply(function(start, end, scale) length(which(!is.na(data[start:end]))) / scale,
                      start=start_ind, end=end_ind, scale=scale)
  
  
  #Set missing cells back to missing
  if (all(is.na(data))) { frequency <- rep(NA, length(frequency)) }
  
  #Return frequency
  return(frequency)
  
}
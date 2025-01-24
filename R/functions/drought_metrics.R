drought_metrics <- function(ts, scale, perc) {
  
  library(zoo)
  
  ts_roll <- rollmean(ts, scale, na.pad=TRUE) 
  
  
  #Calculate climatological values for each month
  climatology <- sapply(1:12, function(x) mean(ts_roll[seq(x, by=12, length.out=length(ts_roll)/12)], na.rm=TRUE))
  
  #Repeat to make same length as time series
  climatology <- rep(climatology, length.out=length(ts))
  
  #Calculate threshold for each month
  threshold <- sapply(1:12, function(x) quantile(ts_roll[seq(x, by=12, length.out=length(ts_roll)/12)], 
                                                 probs=perc, na.rm=TRUE))
  
  ### Indetify months below threshold ###
  drought_months <- ts_roll <= rep(threshold, length.out=length(ts_roll))
  
  
  ### Find drought durations ###
  tmp <- rle(drought_months)
  
  consec <- with(tmp, rep(lengths * values, lengths))
  
  nonzero <- rle(consec != 0) # Logical vector indicating non-zero values
  
  # Extract lengths and positions of consecutive non-zero runs. These are the
  #durations (need to NA omit as these get propagated through)
  duration <- na.omit(nonzero$lengths[nonzero$values == TRUE])
  
  ### Intensity ###
  
  # Find the start indices of all drought events
  run_ends <- cumsum(nonzero$lengths)  # Cumulative end positions of runs
  run_starts <- c(1, head(run_ends, -1) + 1)  # Start positions of all runs
  
  # Extract start indices corresponding to non-zero runs
  non_zero_starts <- na.omit(run_starts[nonzero$values == TRUE])
  
  
  inds <- lapply(1:length(non_zero_starts), function(x) non_zero_starts[x]:(non_zero_starts[x]+duration[x]-1))
  
  intensity <- sapply(inds, function(x) abs(( mean(ts[x]) / mean(climatology[x]) -1)) * 100)
  
  if(any(intensity > 100)) stop("Error calculating intensity")
  
  
  return(list(duration=mean(duration), intensity=mean(intensity)))
  
  
}
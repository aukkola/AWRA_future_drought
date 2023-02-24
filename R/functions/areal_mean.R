areal_mean <- function(hdata) {
  
  out <- vector(length=nlayers(hdata))
  
  for (i in 1:nlayers(hdata)) {
    
    harea=area(hdata[[i]],na.rm=TRUE)
    tarea=sum(values(harea),na.rm=TRUE)
    
    vals <- values(hdata[[i]]*harea)
    
    if (all(is.na(vals))) {
      out[i] <- NA
    } else {
      tot <- sum(vals,na.rm=TRUE)
      
      out[i] <- tot/tarea
      
    }

  }
  
  return(out)
}




areal_mean_terra <- function(data) {
  
  global(data, "mean", na.rm=TRUE, weights=cellSize(data))
  #(cellsize only returns the area for cells that are not NA
  #so don't need to worry about NA-cells here)
}











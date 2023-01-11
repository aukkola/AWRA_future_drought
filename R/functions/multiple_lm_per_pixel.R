multiple_lm_per_pixel <- function(data, npred=2, weights, interaction=FALSE)
{
  
  #Split data to y and predictors (y assumed to come first)
  
  n <- length(data) / (npred + 1)
  
  start_inds <- seq(1, by=n, length.out = npred+1)
  
  
  
  #Dependent data (y)
  y <- data[start_inds[1] : n]
  
  #Independent data (currently set up for 2)
  x1 <- data[start_inds[2] : (n*2)] 
  
  x2 <- data[start_inds[3] : (n*3)] 
  
  
  #Normalise weights
  weigths <- weights / sum(weights)
  
  
  
  #If not data missing
  if (all(!is.na(c(y, x1, x2)))) {
    
    if (interaction) {
      
      stop("not yet implemented")
      
      
    } else {
      
      #Weighted linear regression
      lm <- summary(lm(y ~ x1 + x2, weights=weights))
      
      
      #Collate coefficients and their p-values
      outs <- c(x1=lm$coefficients[2], x2=lm$coefficients[3], 
                p_x1=lm$coefficients[11], p_x2=lm$coefficients[12],
                rsq=lm$r.squared)
      
      
    }
    
    
    return(outs)
    

  } else {
    
    return(rep(NA, 5))
  }
  
 

} #function
  
  
  
  
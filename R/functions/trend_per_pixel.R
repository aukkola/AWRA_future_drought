# Functions --------------------------------------------


#Define function to calculate slope and p-value
trend_per_pixel <- function(data_vec){
  
  len <- length(data_vec)
  
  #Get x and y data
  x_data <- 1:length(data_vec)
  y_data <- data_vec
  
  #All missing, return NA
  if(all(is.na(x_data)) | all(is.na(y_data))) {
    return(c(NA, NA, NA))
    
    
    #Data available, perform linear regression
  } else {
    
    lm <- lm(y_data ~ x_data)
    
    #Get slope, intercept and p-value
    coef <- lm$coefficients[2]
    intercept <- lm$coefficients[1]
    pval <- summary(lm)$coefficients[8]
    
    return(c(coef, intercept, pval))
    
  }
}



#Define function to calculate slope and p-value
#using Mann-Kendall and Sen slope
trend_per_pixel_MannKen <- function(data_vec){
  
  library(trend)

  #All missing, return NA
  if (all(is.na(data_vec))) {
    return(c(NA,NA))
    
    
    #Data available, perform linear regression
  } else {
    
    lm <- sens.slope(data_vec)
    
    #Get slope and p-value
    coef <- lm$estimates
    
    pval <- lm$p.value
    
    return(c(coef,pval))
    
  }
}





  
  
  
  
  

#From package Hmisc, unable to install on raijin
#Original function for variance but taking the square root here to get standard deviation

wtd_stdev <- function(x, weights) { 
  
  stdev <- sqrt(wtd_var(x, weights=weights))
    
  return(stdev)
}


#From Hmisc:
wtd_var <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE, method = c("unbiased", 
                                                                      "ML")) 
{
  method <- match.arg(method)
  if (!length(weights)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    return(var(x))
  }
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  if (normwt) 
    weights <- weights * length(x)/sum(weights)
  if (normwt || method == "ML") 
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))
  sw <- sum(weights)
  if (sw <= 1) 
    warning("only one effective observation; variance estimate undefined")
  xbar <- sum(weights * x)/sw
  sum(weights * ((x - xbar)^2))/(sw - 1)
}
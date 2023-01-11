#Weighted correlation

#function to calculate weighted correlation per raster pixel
wtd_cor_per_pixel <- function(data, weights) {
  
  #Need to split input to x and y
  x <- data[1:(length(data)/2)]
  y <- data[(length(data)/2 +1) : length(data)]
  
  
  #Calculate if all models available
  if(all(!is.na(x)) & all(!is.na(y))) {
    
    #Calculate weighted correlation
    corr <- wtd_cor(x, y, weight=weights)
    
    #Return correlation coefficient and p-value
    return(c(corr=corr[1], pval=corr[4]))
    
  } else {
    return(c(corr=NA, pval=NA))
  }
  
}


#----------------

#These are from package weights

wtd_cor <- function (x, y = NULL, weight = NULL, mean1 = TRUE, collapse = TRUE, 
          bootse = FALSE, bootp = FALSE, bootn = 1000) 
{
  x <- as.matrix(x)
  xnm <- colnames(x)
  if (is.null(weight)) {
    weight <- rep(1, dim(x)[1])
  }
  if (bootse == FALSE & bootp == TRUE) 
    warning("bootp can only be used with bootstrapped standard errors")
  if (mean1 == TRUE) 
    weight <- weight/mean(weight, na.rm = TRUE)
  if (is.null(y)) {
    y <- x
  }
  y <- as.matrix(y)
  ynm <- colnames(y)
  if (is.null(xnm)) 
    xnm <- "X"
  if (is.null(ynm)) 
    ynm <- "Y"
  if (dim(x)[1] != dim(y)[1]) 
    stop("Cannot Correlate Variables of Different Lengths")
  if (bootse == FALSE) {
    materset <- lapply(as.data.frame(x), function(x) lapply(as.data.frame(y), 
                                                            function(y) onecor.wtd(x, y, weight)))
    est <- sapply(materset, function(q) sapply(q, function(g) g[1]))
    se <- sapply(materset, function(q) sapply(q, function(g) g[2]))
    tval <- sapply(materset, function(q) sapply(q, function(g) g[3]))
    pval <- sapply(materset, function(q) sapply(q, function(g) g[4]))
    out <- list(correlation = est, std.err = se, t.value = tval, 
                p.value = pval)
  }
  if (bootse == TRUE) {
    est <- as.matrix(wtd.cors(x, y, weight))
    samps <- lapply(1:bootn, function(g) sample(1:dim(x)[1], 
                                                round(sum(weight, na.rm = TRUE), 0), replace = TRUE, 
                                                prob = weight))
    corset2 <- lapply(samps, function(q) as.matrix(cor(x[q, 
                                                         ], y[q, ], use = "pairwise.complete.obs")))
    eachcor <- lapply(1:dim(est)[1], function(a) sapply(1:dim(est)[2], 
                                                        function(b) unlist(sapply(corset2, function(g) g[a, 
                                                                                                         b]))))
    est2 <- sapply(eachcor, function(a) colMeans(a))
    se <- sapply(eachcor, function(a) sqrt(apply(a, 2, var)))
    tval <- est2/se
    pval <- pchisq(tval^2, 1, lower.tail = FALSE)
    if (bootp == TRUE) 
      pval <- sapply(eachcor, function(a) apply(a, 2, function(x) 2 * 
                                                  min(c(sum(x > 0 & !is.na(x))/sum(!is.na(x)), 
                                                        sum(x < 0 & !is.na(x))/sum(!is.na(x))))))
    if (length(ynm) > 1 & length(xnm) > 1) {
      colnames(est2) <- colnames(se) <- colnames(tval) <- colnames(pval) <- xnm
      rownames(est2) <- rownames(se) <- rownames(tval) <- rownames(pval) <- ynm
    }
    out <- list(correlation = t(est), bootcor = est2, std.err = se, 
                t.value = tval, p.value = pval)
  }
  if (is.vector(est) & collapse == TRUE || (1 %in% dim(est)) & 
      collapse == TRUE) {
    outpre <- out
    if (bootse == FALSE) 
      out <- matrix(unlist(out), ncol = 4, byrow = FALSE)
    if (bootse == TRUE) 
      out <- matrix(unlist(out), ncol = 5, byrow = FALSE)
    nom <- xnm
    if (length(xnm) == 1) 
      nom <- ynm
    rownames(out) <- nom
    colnames(out) <- names(outpre)
  }
  out
}

#-----------------------------

onecor.wtd <- function(x, y, weight=NULL){
  if(sum((!is.na(x))*(!is.na(y)))>0){
    if(is.null(weight)){
      weight <- rep(1, length(x))
    }
    use <- !is.na(y) & !is.na(x)
    x <- x[use]
    y <- y[use]
    weight <- weight[use]
    #r1 <- lm(stdz(y, weight=weight)~stdz(x, weight=weight), weight=weight)
    #corcoef <- coef(summary(r1))[2,]
    corcoef <- coef(summary(lm(stdz(y, weight=weight)~stdz(x, weight=weight), weights=weight)))[2,]
  }
  else
    corcoef <- rep(NA, 4)
  names(corcoef) <- rep("", 4)
  corcoef
}

#-----------------------------

stdz <- function (x, weight = NULL) 
{
  if (is.null(weight)) {
    weight <- rep(1, length(x))
  }
  x <- x - wtd.mean(x, weight, na.rm = TRUE)
  x <- x/sqrt(wtd.var(x, weight, na.rm = TRUE))
  x
}

#-----------------------------

wtd.mean <- function (x, weights = NULL, normwt = "ignored", na.rm = TRUE) 
{
  if (!length(weights)) 
    return(mean(x, na.rm = na.rm))
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  sum(weights * x)/sum(weights)
}

#-----------------------------

wtd.var <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE, method = c("unbiased", 
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
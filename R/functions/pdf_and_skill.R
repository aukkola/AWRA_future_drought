perkins_skill_score <- function(model, obs, nbins=50) {

  # #First calculate equal breaks (R hist function doesn't let you
  # #match the number of bins and their sizes nicely without doing this manually)
  # range <- range(c(model, obs), na.rm=TRUE)
  # 
  # bins <- seq(range[1], range[2], length.out=nbins)
  
  #Calculate skill score
  skill <- sum(sapply(1:length(model), function(x) min(model[x], obs[x])), na.rm=TRUE)
  
  if(skill < 0 | skill >1) stop("wrong Perkins skill value determined")
  
  return(skill)
}

perkins_skill_score_normalised <- function(model, obs) {
  
  
  #First normalise probability density values to make sure skill score is 
  #between 0-1
  model <- model/sum(model)
  obs   <- obs / sum(obs)
  
  
  #Calculate skill score
  skill <- sum(sapply(1:length(model), function(x) min(model[x], obs[x])), na.rm=TRUE)
  
  if(skill < 0 | skill >1) stop("wrong Perkins skill value determined")
  
  return(skill)
}


#Function to calculate probability densities
hist_xy <- function(data, n=50) {
  
  hist <- hist(data, breaks=n, plot=FALSE)
  #   offset <- hist$breaks[2] - hist$breaks[1]
  #  density <- hist$density*offset #normalise density to sum up to 1
  return(hist$counts / sum(hist$counts))
}


# #Function to calculate probability densities
# hist_xy <- function(data, n=50) {
#   
#   hist <- hist(data, breaks=n, plot=FALSE)
#   #   offset <- hist$breaks[2] - hist$breaks[1]
#   #  density <- hist$density*offset #normalise density to sum up to 1
#   return(list(x=hist$breaks, y=hist$density))
# }


# Metric_NME.R
#
# Gab Abramowitz UNSW 2020 (gabsun at gmail dot com)
#grabbed from here: https://gitlab.com/modelevaluation/me.org-r-library/-/blob/master/pals/R/Metric_NME.R?ref_type=heads
NME = function(model,obs){
  if(all(is.na(model)) | all(is.na(obs))){
    NME = NA
  }else{
    NME = sum(abs(as.vector(model)-as.vector(obs)),na.rm=TRUE)/sum(abs( mean(obs,na.rm=TRUE) -obs), na.rm=TRUE)
  }
  return(NME)
}


#Calculate fraction of land area in each histogram bin
#Function from GCB paper
freqs <- function(x, breaks, area) {
  
  area_tot <- sum(area, na.rm=TRUE)
  
  density <- vector()
  
  for (b in 1:(length(breaks)-1)) {
    
    ind <- which(x >= breaks[b] & x < breaks[b+1])
    
    density[b] <- sum(area[ind], na.rm=TRUE) / area_tot
  }
  
  # h <- hist(x, breaks=breaks, plot=FALSE)
  # density <- h$counts / sum(h$counts)
  return(density)
}



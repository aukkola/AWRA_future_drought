library(ncdf4)
library(raster)


#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w97/amu561/Steven_CABLE_runs/"


scale=3

#plot difference in historical and future periods
#(using 36 years)
# historical is 1960-2005 but ignoring first 10 years for spinup
#leaving 36yrs of data

all_years <- 1960:2099

hist_yrs <- c(1970, 2005) #start and end 
fut_yrs  <- c(2064, 2099)

#Variables to plot  
variables <- c("pr", "qtot", "sm")

rcps <- c("rcp45", "rcp85")

#Get bias correction methods
bc_methods <- list.files(paste0(path, "/drought_metrics/", scale, "-month/"))


#Initialise
duration <- list()
intensity <- list()


#RCPs
for (r in 1:length(rcps)) {
  
  #Initialise
  duration[[r]] <- list()
  intensity[[r]] <- list()
  
  #Variables
  for (v in 1:length(variables)) {
    
    
    #Save name
    exp_name <- vector()
    
  #BC methods
  for (b in 1:length(bc_methods)) {
    
    #Get GCMs
    gcms <- list.files(paste0(path, "/drought_metrics/", scale, "-month/", bc_methods[b]))
     
    #GCMs
    for (g in 1:length(gcms)) {
      
         
        #Initialise
        if (b==1 & g==1) {
          duration[[r]][[v]]  <- brick()
          intensity[[r]][[v]] <- brick()
          
        }
        

        ################
        ### Get data ###
        ################
        
        #Find data file
        data_file <- list.files(paste0(path, "/drought_metrics/", scale, "-month/", 
                                       bc_methods[b], "/", gcms[g]), 
                                pattern=paste(variables[v], rcps[r], scale, sep="_"),
                                full.names=TRUE)
        
        if(length(data_file) ==0 ) {
          
          print(paste0("couldn't find a file, r:", r, ", b:", b, ", g:", g,
                       ", v:", v))
          next
        }
        
        #Save name
        exp_name <- append(exp_name, paste0(bc_methods[b], "_", gcms[g]))
        
        
        #Get metrics data
        diff_duration  <- raster(data_file, varname="duration")
        diff_intensity <- raster(data_file, varname="intensity")
      
        
        
        #Save to list
        duration[[r]][[v]]  <- addLayer(duration[[r]][[v]], diff_duration)
        intensity[[r]][[v]] <- addLayer(intensity[[r]][[v]], diff_intensity)
        
        
      
      
    } #GCMs
    
  } #BC
    
    names(duration[[r]][[v]]) <- exp_name
    
    
  } #variables
  names(duration[[r]]) <- variables
  
} #RCPs


names(duration)  <- rcps
names(intensity) <- rcps


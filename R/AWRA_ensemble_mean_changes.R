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
  
  
  #BC methods
  for (b in 1:length(bc_methods)) {
    
    #Get GCMs
    gcms <- list.files(paste0(path, "/drought_metrics/", scale, "-month/", bc_methods[b]))
    
    
    #GCMs
    for (g in 1:length(gcms)) {
      
      #Variables
      for (v in 1:length(variables)) {
        
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
        
        #Get metrics data
        temp_duration  <- brick(data_file, varname="duration")
        temp_intensity <- brick(data_file, varname="intensity")
        
        
        #Check that both data have the correct number of layers
        if (nlayers(temp_duration) != (length(all_years)*12) |
            nlayers(temp_intensity) != (length(all_years)*12)) {
          stop("something funny going on")
        }
        
        
        #Calculate historical and future mean
        ind_hist <- c(which(all_years == hist_yrs[1])*12-11, #Jan of first historical yr
                      which(all_years == hist_yrs[2])*12)
        
        ind_fut <- c(which(all_years == fut_yrs[1])*12-11, #Jan of first future yr
                     which(all_years == fut_yrs[2])*12)
        
        
        hist_duration <- mean(temp_duration[[ind_hist[1]:ind_hist[2]]], na.rm=TRUE)
        fut_duration  <- mean(temp_duration[[ind_fut[1]:ind_fut[2]]], na.rm=TRUE)
        
        hist_intensity <- mean(temp_intensity[[ind_hist[1]:ind_hist[2]]], na.rm=TRUE)
        fut_intensity  <- mean(temp_intensity[[ind_fut[1]:ind_fut[2]]], na.rm=TRUE)
        
        
        #Calculate difference future - historical
        diff_duration  <- fut_duration - hist_duration
        diff_intensity <- fut_intensity - hist_intensity
        

        #Save to list
        duration[[r]][[v]]  <- addLayer(duration[[r]][[v]], diff_duration)
        intensity[[r]][[v]] <- addLayer(intensity[[r]][[v]], diff_intensity)
        
        
      } #variables
      
    } #GCMs
    
  } #BC
  
} #RCPs


#Save RDS files so don't need to run this processing again to plot

out_path <- paste0(path, "/Mean_drought_metrics")

dir.create(out_path)

#Duration
saveRDS(duration, paste0(out_path, "/AWRA_future_difference_in_duration_scale_", scale,
                         "_", hist_yrs[1], "-", hist_yrs[2], "_vs_", fut_yrs[1],
                         "-", fut_yrs[2], "_all_variables_and_RCPs.nc"))

#Duration
saveRDS(intensity, paste0(out_path, "/AWRA_future_difference_in_intensity_scale_", scale,
                          "_", hist_yrs[1], "-", hist_yrs[2], "_vs_", fut_yrs[1],
                          "-", fut_yrs[2], "_all_variables_and_RCPs.nc"))





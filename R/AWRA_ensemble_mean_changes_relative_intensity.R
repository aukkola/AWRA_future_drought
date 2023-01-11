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

#Output path
out_path <- paste0(path, "/Mean_drought_metrics/")



#RCPs
for (r in 1:length(rcps)) {
  
  
  #Variables
  for (v in 1:length(variables)) {

    out_path_intensity <- paste(out_path, rcps[r], variables[v], "rel_intensity", sep="/")
    dir.create(out_path_intensity, recursive=TRUE)
    
    #Save name
    #exp_name <- vector()
    
    
    #BC methods
    for (b in 1:length(bc_methods)) {
      
      #Get GCMs
      gcms <- list.files(paste0(path, "/drought_metrics/", scale, "-month/", bc_methods[b]))
      
      
      #GCMs
      for (g in 1:length(gcms)) {
        
    
        ### Output file names ###
        
 
        out_file_intensity <- paste0(out_path_intensity, "/AWRA_future_difference_in_rel_intensity_scale_", 
                                     scale, "_", hist_yrs[1], "-", hist_yrs[2], "_vs_", fut_yrs[1],
                                     "-", fut_yrs[2], "_", variables[v], "_", rcps[r], "_",
                                     bc_methods[b], "_", gcms[g], ".nc")
        
        
        #If these exist already, skip to speed up code
        
        if (file.exists(out_file_intensity)) {
          print ("skipping")
          next
        }
        
        
        ################
        ### Get data ###
        ################
        
        #Find data file
        data_file <- list.files(paste0(path, "/drought_metrics/", scale, "-month/", 
                                       bc_methods[b], "/", gcms[g]), 
                                pattern=paste(variables[v], rcps[r], scale, sep="_"),
                                full.names=TRUE)
        
        #A few files missing currently, skip if not found
        if(length(data_file) ==0 ) {
          print(paste0("couldn't find a file, r:", r, ", b:", b, ", g:", g,
                       ", v:", v))
          next
        }
        
        
        #Get metrics data
        temp_intensity <- brick(data_file, varname="rel_intensity")
        
        
        #Check that both data have the correct number of layers
        if (nlayers(temp_intensity) != (length(all_years)*12)) {
          stop(paste0("something funny going on, r:" 
                      , r, ", b:", b, ", g:", g,", v:", v))
        }
        
        
        #Save name
        #exp_name <- append(exp_name, paste0(bc_methods[b], "_", gcms[g]))
        
        
        #Calculate historical and future mean
        ind_hist <- c(which(all_years == hist_yrs[1])*12-11, #Jan of first historical yr
                      which(all_years == hist_yrs[2])*12)
        
        ind_fut <- c(which(all_years == fut_yrs[1])*12-11, #Jan of first future yr
                     which(all_years == fut_yrs[2])*12)
        
        
        hist_intensity <- mean(temp_intensity[[ind_hist[1]:ind_hist[2]]], na.rm=TRUE)
        fut_intensity  <- mean(temp_intensity[[ind_fut[1]:ind_fut[2]]], na.rm=TRUE)
        
        
        #Calculate difference future - historical
        diff_intensity <- fut_intensity - hist_intensity
        
        
        
        #Write netcdf
        writeRaster(diff_intensity, out_file_intensity, varname="intensity", overwrite=TRUE)
        
        
        
        
      } #GCMs
      
    } #BC
    
    #names(duration[[r]][[v]]) <- exp_name
    #names(intensity[[r]][[v]]) <- exp_name
    
    
  } #variables
  
  #names(duration[[r]])  <- variables
  #names(intensity[[r]]) <- variables
  
} #RCPs


#names(duration)  <- rcps
#names(intensity) <- rcps



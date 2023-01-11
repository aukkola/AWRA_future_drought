calculate_time_mean <- function(models, outdir_exp,  exp, proj, metrics, extent_raster)
{
  
  #Not yet implemented for frequency, stop
  if (metrics == "frequency") { stop("Not implemented for frequency") }
  
  
  #Set period for calculating mean
  if (exp %in% c("historical", "observed")) {
    out_period <- "1950_2014"
  } else {
    out_period <- "2051_2099"
    
  }
  
  
  #Find first ensemble member
  files <- list.files(models, full.names=TRUE)
  
  
  
  for (f in 1:length(files)) {
    
    
    print(paste0("model:", models))
    
    
    #Set output file name and check if it already exists
    outfile <- paste0(outdir_exp, "/Mean_changes_", out_period, "_", basename(files[f]))
    
    
    
    if(metrics != "frequency") { if (file.exists(outfile)) { next } }
    
    
    if (metrics == "frequency") {
      data <- brick(files[f], varname="duration", stopIfNotEqualSpaced=FALSE)
      
    } else {
      data <- brick(files[f], varname=metrics, stopIfNotEqualSpaced=FALSE)
    }
    
    
    #Get time period
    data_yrs <- as.numeric(format(as.Date(data@z$Date), "%Y"))
    
    
    
    #REMOVE ONCE FIXED IN PYTHON CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #Years are written incorrectly to some python-produced files
    #The data should be fine though so correct years here
    if (exp == "historical"  & proj == "CMIP5" & data_yrs[length(data_yrs)] < 2014) {
      
      ind <- tail(which(data_yrs==2005), n=1) + 1
      
      data_yrs[ind:length(data_yrs)] <- rep(2006:2014, each=12)
      
    }
    
    
    
    #Extract time periods
    
    #If historical or obs, set years to 1950-2014
    if (exp %in% c("historical", "observed")) {
      
      #Need to add exception to CESM1-WACCM, only reports outputs from 1955
      if (grepl("CESM1-WACCM", models)) {
        start_yr <- 1955
      } else {
        start_yr <- 1950
      }
      
      if (metrics == "frequency") {
        
        #Express as events per 120 months (i.e. 10 years)
        subset_data <- calc(data[[which(data_yrs == start_yr)[1] : tail(which(data_yrs == 2014), n=1)]], 
                            fun=calc_freq)
        
      } else {
        
        subset_data <- mean(data[[which(data_yrs == start_yr)[1] : 
                                    tail(which(data_yrs == 2014), n=1)]], na.rm=TRUE)
      }
      
      
      
      #And two periods for future (2015-2050, 2051-2099)
    } else {
      
      
      if (metrics == "frequency") {
        
        #Express as events per 120 months (i.e. 10 years)
        
        #2015-2050
        subset_data <- calc(data[[which(data_yrs == 2015)[1] : tail(which(data_yrs == 2050), n=1)]],
                            fun=calc_freq) #ction(x) length(which(!is.na(x))) / length(x))
        
        #2051-2099
        subset_data <- addLayer(subset_data, calc(data[[which(data_yrs == 2051)[1] : tail(which(data_yrs == 2099), n=1)]],
                                                  fun=calc_freq)) #ction(x) length(which(!is.na(x)))) / length(x))
        
        
      } else {
        
        #2015-2050
        subset_data <- mean(data[[which(data_yrs == 2015)[1] : 
                                    tail(which(data_yrs == 2050), n=1)]], na.rm=TRUE)
        
        #2051-2099 (use 2099 because some CMIP5 models finish then)
        subset_data <- addLayer(subset_data, mean(data[[which(data_yrs == 2051)[1] : 
                                                          tail(which(data_yrs == 2099), n=1)]], na.rm=TRUE))
        
      }
      
      
      
    }
    
    
    #Need to mask frequency
    if (metrics == "frequency") { 
      subset_data <- mask(subset_data, mean(data, na.rm=TRUE))
    }
    
    
    #Resample data and calculate mean
    resample_data <- resample(subset_data, extent_raster)
    
    
    #Write output
    writeRaster(resample_data, outfile, varname=metrics, longname=paste("Trend in", metrics, out_period),
                overwrite=TRUE)
    
  } #files
  
  
  
  
  
}   
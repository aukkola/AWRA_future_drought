

#Daily flt files 
write_daily_flt <- function(data, variable, outdir) {
  
  #Get dates
  dates <- getZ(data)
  
  #Write daily flt files
  for (n in 1:nlayers(data)) {
    
    day_stamp <- format(dates[n], format="%Y%m%d")
    
    #Following Mengyuan's naming convention
    outfile <- paste0(outdir, "/", variable, "_", day_stamp, "_",
                      day_stamp, ".flt")
    
    #Write output
    writeRaster(data[[n]], outfile, overwrite=TRUE)
    
  }
  
}



#Daily LWdown flt files 
write_daily_flt_lwdown <- function(tmean, tmin, swdown, variable, outdir) {
  
  #Get dates
  dates <- getZ(tmean)
  
  #Day of year
  doy <- as.numeric(strftime(dates, format = "%j"))
  
  #Latitude
  lat <- raster::init(tmean, 'y')
  
  
  #Write daily flt files
  for (n in 1:length(dates)) {
    
    
    ### Calculate LWdown for the day ###
    
    #Stack inputs for calc function
    ins <- brick(tmean[[n]], tmin[[n]], swdown[[n]], lat)
    
    #Calculate LWdown
    lwdown <- calc(ins, function(x) calculate_lwdown(tmean_K=x[1], 
                                                     tmin_K=x[2], 
                                                     swdown_Wm2=x[3], 
                                                     latitude_deg=x[4], 
                                                     DOY=doy[n]))
    
    ### Write output ###
    
    day_stamp <- format(dates[n], format="%Y%m%d")
    
    #Following Mengyuan's naming convention
    outfile <- paste0(outdir, "/", variable, "_", day_stamp, "_",
                      day_stamp, ".flt")
    
    #Write output
    writeRaster(lwdown, outfile, overwrite=TRUE)
    
  }
  
}  




#Daily VPD flt files 
write_daily_flt_vpd <- function(tmean, tmin, variable, outdir) {
  
  
  #Get dates
  dates <- getZ(tmean)
  
  
  #Write daily flt files
  for (n in 1:length(dates)) {
    
    
    ### Calculate VPD for the day ###
    
    #Stack inputs for calc function
    ins <- brick(tmean[[n]], tmin[[n]])
    
    
    #Calculate VPD using constant air pressure
    vpd <- calc(ins, function(x) calculate_vpd(tmean=x[1], 
                                                 tmin=x[2]))
    
    
    ### Write output ###
    
    day_stamp <- format(dates[n], format="%Y%m%d")
    
    #Following Mengyuan's naming convention
    outfile <- paste0(outdir, "/", variable, "_", day_stamp, "_",
                      day_stamp, ".flt")
    
    #Write output
    writeRaster(vpd, outfile, overwrite=TRUE)
    
    
  }
  
}  



#Daily flt files 
write_daily_flt <- function(data, variable, outdir) {
  
  #Get dates
  dates <- getZ(data)
  
  #Write daily flt files
  for (n in 1:nlayers(data)) {
    
    #Skip years between 2006-2053 as only running CABLE for 2054-2099
    yr <- as.numeric(format(dates[n], format="%Y"))
    
    #if (yr >2005 & yr < 2054) next
    
    
    day_stamp <- format(dates[n], format="%Y%m%d")
    
    #Following Mengyuan's naming convention
    outfile <- paste0(outdir, "/", variable, "_", day_stamp, "_",
                      day_stamp)
    
    #Add this to avoid duplicate processing
    if (file.exists(outfile)) next
    
    #Write output
    writeRaster(data[[n]], paste0(outfile, ".flt"), overwrite=TRUE)
    
    
    #Remove extra files, not needed by weather generator
    unlink(paste0(outfile, ".flt.aux.xml"))
    unlink(paste0(outfile, ".hdr"))
    unlink(paste0(outfile, ".prj"))
    
    
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
    
    #Skip years between 2006-2053 as only running CABLE for 2054-2099
    yr <- as.numeric(format(dates[n], format="%Y"))
    
    #if (yr >2005 & yr < 2054) next
    
    
    
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
                      day_stamp)
    
    #Add this to avoid duplicate processing
    if (file.exists(outfile)) next
    
    #Write output
    writeRaster(lwdown, paste0(outfile, ".flt"), overwrite=TRUE)
    
    
    #Remove extra files, not needed by weather generator
    unlink(paste0(outfile, ".flt.aux.xml"))
    unlink(paste0(outfile, ".hdr"))
    unlink(paste0(outfile, ".prj"))
    
    
  }
  
}  




#Daily VPD flt files 
write_daily_flt_vpd <- function(tmax, tmin, variable, outdir) {
  
  
  #Get dates
  dates <- getZ(tmax)
  
  
  #Write daily flt files
  for (n in 1:length(dates)) {
    
    
    #Skip years between 2006-2053 as only running CABLE for 2054-2099
    yr <- as.numeric(format(dates[n], format="%Y"))
    
   # if (yr >2005 & yr < 2054) next
    
    
    ### Calculate VPD for the day ###
    
    #Stack inputs for calc function
    ins <- brick(tmax[[n]], tmin[[n]])
    
    
    #Calculate VPD using constant air pressure
    vpd <- calc(ins, function(x) calculate_vpd(tmax=x[1], tmin=x[2]))
    
    
    ### Write output ###
    
    day_stamp <- format(dates[n], format="%Y%m%d")
    
    #Following Mengyuan's naming convention
    outfile <- paste0(outdir, "/", variable, "_", day_stamp, "_",
                      day_stamp)
    
    #Add this to avoid duplicate processing
    if (file.exists(outfile)) next
    
    #Write output
    writeRaster(vpd, paste0(outfile, ".flt"), overwrite=TRUE)
    
    
    #Remove extra files, not needed by weather generator
    unlink(paste0(outfile, ".flt.aux.xml"))
    unlink(paste0(outfile, ".hdr"))
    unlink(paste0(outfile, ".prj"))
    
  }
  
}  

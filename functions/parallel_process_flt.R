
parallel_process_flt <- function(variable, outdir_all, awra_path, model, 
                                 experiment, bc_method) {
  
  
  #Output directory name
  outdir<- paste0(outdir_all, "/", variable)
  
  dir.create(outdir, recursive=TRUE)
  
  
  ### LWdown ###
  
  #(not used, calculated using weather generated data)
  
  if(variable == "lwdown") {
    
    #Get Tmean
    tmean <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                    bc_method, "latest/day/tas/", sep="/"),
                              full.names=TRUE), varname="tas")
    
    #Get Tmin
    tmin <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                   bc_method, "latest/day/tasmin/", sep="/"),
                             full.names=TRUE), varname="tasmin")
    
    #Get SWdown
    swdown <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                     bc_method, "latest/day/rsds/", sep="/"),
                               full.names=TRUE), varname="rsds")
    
    
    
    #Write daily flts
    write_daily_flt_lwdown(tmean=tmean, tmin=tmin, swdown=swdown, 
                           variable=variable, outdir=outdir)
    
    
    
    
  ### VPD ###
    
  } else if (variable %in% c("vpd_09", "vpd_15")) {
    
    
    if (variable == "vpd15") {
      
      #Get Tmax for 3pm VPD
      tmax <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                      bc_method, "latest/day/tasmax/", sep="/"),
                                full.names=TRUE), varname="tasmax")
      
      write_daily_flt_vpd(tair=tmax, variable=variable,
                          outdir=outdir)
      
      
      
    } else if (variable == "vpd09") {
      
      #Get Tmin for 9am VPD
      tmin <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                     bc_method, "latest/day/tasmin/", sep="/"),
                               full.names=TRUE), varname="tasmin")
      
      write_daily_flt_vpd(tair=tmin, variable=variable,
                          outdir=outdir)
      
      
    }
    
     
  
  ### All other variables ###
      
  } else {
    
    #Find data file  
    data_file <- list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                  bc_method, "latest/day", variable, sep="/"),
                            full.names=TRUE)
    
    #Read data
    data <- brick(data_file, varname=variable)
    
    
    #Write daily flts
    write_daily_flt(data=data, variable=variable, 
                    outdir=outdir)
    
  }
  
}




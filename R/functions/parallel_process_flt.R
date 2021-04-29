
parallel_process_flt <- function(variable, outdir_all, awra_path, model, 
                                 experiment, bc_method) {
  
  
  #Output directory name
  outdir<- paste0(outdir_all, "/", variable)
  
  dir.create(outdir, recursive=TRUE)
  
  print(paste0("variable: ", variable))
  
  
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
    
    #Use Tmax to calculate 3pm saturated vapour pressure and
    #Tmean to calculate 9am saturated vapour pressure
    #Tmin always used to calculate actual vapour pressure
    #otherwise difficult to create a diurnal cycle
    
    if (variable == "vpd_15") {
      
      #Get Tmax for saturated VP
      tmax <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                      bc_method, "latest/day/tasmax/", sep="/"),
                                full.names=TRUE), varname="tasmax")
      
      #Get Tmin for actual VP
      tmin <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                     bc_method, "latest/day/tasmin/", sep="/"),
                               full.names=TRUE), varname="tasmin")
      
      
      
      write_daily_flt_vpd(tmax=tmax, tmin=tmin, variable=variable,
                          outdir=outdir)
      
      
      
    } else if (variable == "vpd_09") {
      
      #Get Tmean for saturated VP
      tmean <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                     bc_method, "latest/day/tas/", sep="/"),
                               full.names=TRUE), varname="tas")
      
      #Get Tmin for actual VP
      tmin <- brick(list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                     bc_method, "latest/day/tasmin/", sep="/"),
                               full.names=TRUE), varname="tasmin")
      
      
      write_daily_flt_vpd(tmax=tmean, tmin=tmin, variable=variable,
                          outdir=outdir)
      
      
    } else {
      stop("VPD variable name not recognised")
    }
    
     
  
  ### All other variables ###
      
  } else {
    
    #Find data file  
    data_file <- list.files(paste(awra_path, model, experiment, "r1i1p1", 
                                  bc_method, "latest/day", variable, sep="/"),
                            full.names=TRUE)
    
    #Read data
    data <- brick(data_file, varname=variable)
    
    #Not needed, doing this in fortran code
    # #Convert temperature from K to C
    # if (variable %in% c("tas", "tasmax", "tasmin")) {
    #   
    #   data <- data - 273.15
    #   
    # }
    # 
    # #Convert rainfall from mm/s to mm/day
    # if (variable == "pr") {
    #   
    #   data <- data * 86400
    #   
    # }
    # 
    # #Convert SWdown from W m-2 to MJ/day
    # if(variable == "rsds") {
    #   
    #   data <- data * 86400 / 10**6
    #   
    # }
    # 
    
    #Write daily flts
    write_daily_flt(data=data, variable=variable, 
                    outdir=outdir)
    
  }
  
}




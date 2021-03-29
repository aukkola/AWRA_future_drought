
parallel_process_flt <- function(variable, outdir_all) {
  
  
  ### LWdown ###
  
  if(variable == "lwdown") {
    
    #Get Tmean
    tmean <- brick(list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                    bc_methods[b], "latest/day/tas/", sep="/"),
                              full.names=TRUE), varname="tas")
    
    #Get Tmin
    tmin <- brick(list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                   bc_methods[b], "latest/day/tasmin/", sep="/"),
                             full.names=TRUE), varname="tasmin")
    
    #Get SWdown
    swdown <- brick(list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                     bc_methods[b], "latest/day/rsds/", sep="/"),
                               full.names=TRUE), varname="rsds")
    
    
    
    #Write daily flts
    write_daily_flt_lwdown(tmean=tmean, tmin=tmin, swdown=swdown, 
                           variable=variable, outdir=paste0(outdir_all, "/", variable))
    
    
    
    
    ### VPD ###
    
  } else if (variable=="vpd") {
    
    
    #Get Tmean
    tmean <- brick(list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                    bc_methods[b], "latest/day/tas/", sep="/"),
                              full.names=TRUE), varname="tas")
    
    #Get Tmin
    tmin <- brick(list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                   bc_methods[b], "latest/day/tasmin/", sep="/"),
                             full.names=TRUE), varname="tasmin")
    
    
    
    
    write_daily_flt_vpd(tmean=tmean, tmin=tmin, variable=variable,
                        outdir=paste0(outdir_all, "/", variable))
      
      
      
      
      ### All other variables ###
      
  } else {
    
    #Find data file  
    data_file <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                  bc_methods[b], "latest/day", variables[v], sep="/"),
                            full.names=TRUE)
    
    #Read data
    data <- brick(data_file, varname=variables[v])
    
    #Write daily flts
    write_daily_flt(data=data, variable=variable, 
                    outdir=paste0(outdir_all, "/", variable))
    
  }
  
}




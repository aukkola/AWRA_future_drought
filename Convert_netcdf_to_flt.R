## NEED TO LOAD GDAL !!!!!!

library(raster)
library(rgdal)
library(parallel )

#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w35/amu561/Steven_CABLE_runs"

awra_path <- "/g/data/wj02/COMPLIANT/HMINPUT/output/AUS-5/BoM/"








#Use these options to read AWRA inputs
# DOMAIN = AUS-5
# ENSMEMBER = r1i1p1
# VERSION = latest


#Extra variables that need to be calculated
#(weather generator takes VPD (hPa) as input, converts it to Qair)
extra_vars <- c("lwdown", "vpd")




#####################
### Process files ###
#####################


#Initialise cores
cl <- makeCluster(getOption('cl.cores', 8))



### list models ###

models <- list.files(awra_path)

clusterExport(cl, 'models')


#Loop through models
for (m in 1:length(models)) {
  
  #Get experiments
  experiments <- list.files(paste0(awra_path, "/", models[m]))

  clusterExport(cl, 'experiments')
  
  
  #Loop through experiments
  for (e in 1:length(experiments)) {
    
    #Get bias correction methods
    bc_methods <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", sep="/"))
      
    clusterExport(cl, 'bc_methods')
    
    
    #Output directory for all variables
    outdir_all <- paste(path, "CABLE_inputs/Weather_generator_inputs", models[m],
                         experiments[e], bc_methods[b], sep="/")
    
    #Create directory
    dir.create(outdir_all, recursive=TRUE)
      
    
    #Loop through BC methods
    for (b in 1:length(bc_methods)) {
      
      variables <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                              bc_methods[b], "latest/day", sep="/"))

      #Add extra variables
      variables <- append(variables, extra_vars)
      

      clusterExport(cl, list('variables', 'm', 'b', 'e',
                            'parallel_process_flt',
                            'write_daily_flt_lwdown',
                            'write_daily_flt_vpd',
                            'write_daily_flt',
                            'calculate_vpd',
                            'calculate_lwdown'))
      
      
      #Parallel process variables
      parLapply(cl, variables, function(x) parallel_process_flt(variable=x, 
                                                                outdir_all=outdir_all))
      
      
         
    } #bc methods

  } #experiments
  
} #models










# 
# 
# coords <- matrix(c(148.1517, -35.6566), ncol=2)
# 
# #swdown
# v=2
# data_file <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
#                               bc_methods[b], "latest/day", variables[v], sep="/"),
#                         full.names=TRUE)
# 
# data <- brick(data_file, varname=variables[v])
# 
# swdown=as.vector(extract(data, coords))
# 
# 
# #tmean
# v=4
# data_file <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
#                               bc_methods[b], "latest/day", variables[v], sep="/"),
#                         full.names=TRUE)
# 
# data <- brick(data_file, varname=variables[v])
# 
# tmean=as.vector(extract(data, coords))
# 
# 
# #tmin
# v=6
# 
# data_file <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
#                               bc_methods[b], "latest/day", variables[v], sep="/"),
#                         full.names=TRUE)
# 
# data <- brick(data_file, varname=variables[v])
# 
# tmin=as.vector(extract(data, coords))
# 
# 
# 
# 
# 
# lwdown <- calculate_lwdown(tmean_K=tmean[1:365], tmin_K=tmin[1:365], swdown_Wm2=swdown[1:365], 
#                            latitude_deg=coords[1,2], DOY=1:365)
# 






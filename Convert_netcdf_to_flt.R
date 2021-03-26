## NEED TO LOAD GDAL !!!!!!

library(raster)
library(rgdal)


#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w35/amu561/Steven_CABLE_runs"


awra_path <- "/g/data/wj02/COMPLIANT/HMINPUT/output/AUS-5/BoM/"

# DOMAIN = AUS-5
# ENSMEMBER = r1i1p1
# VERSION = latest


### list models ###

models <- list.files(awra_path)


#Loop through models
for (m in 1:length(models)) {
  
  experiments <- list.files(paste0(awra_path, "/", models[m]))

  
  #Loop through experiments
  for (e in 1:length(experiments)) {
    
    bc_methods <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", sep="/"))
      
    
    #Loop through BC methods
    for (b in 1:length(bc_methods)) {
      
      variables <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                              bc_methods[b], "latest/day", sep="/"))
      
      
      for (v in 1:length(variables)) {
        
        data_file <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", 
                                      bc_methods[b], "latest/day", variables[v], sep="/"),
                                full.names=TRUE)
        
        data <- brick(data_file, varname=variables[v])
        
        
        
        writeRaster(data, paste0(path, "/test.flt"))
        
        
      } #variables
      
    } #bc methods

  } #experiments
  
} #models




write_annual_flt <- function(data, variable) {
  
  
  if (variable == "LWdown") {
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
}












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






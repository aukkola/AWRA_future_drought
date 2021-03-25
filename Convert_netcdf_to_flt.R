library(raster)

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
  
  experiments <- list.files(paste0(path, "/", models[m]))

  
  #Loop through experiments
  for (e in 1:length(experiments)) {
    
    bc_methods <- list.files(paste(path, models[m], ensembles[e], "r1i1p1", sep="/"))
      
    
    #Loop through BC methods
    for (b in 1:length(bc_methods)) {
      
      data_file <- list.files(paste(path, models[m], ensembles[e], "r1i1p1", 
                              bc_methods[b], "latest/day", sep="/"), full.names=TRUE)
      
      
      
    }

  }
  
}


  





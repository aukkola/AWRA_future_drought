## NEED TO LOAD GDAL !!!!!!

library(raster)
library(rgdal)
library(parallel )

#clear R environment
rm(list=ls(all=TRUE))

#Command line args
args = commandArgs(trailingOnly=TRUE)

#for testing
# args=c('/g/data/w35/amu561/Steven_CABLE_runs',
#        '/g/data/wj02/COMPLIANT/HMINPUT/output/AUS-5/BoM/',
#        '/scratch/w35/amu561/Steven_CABLE_runs CNRM-CERFACS-CNRM-CM5',
#        'historical',
#        'CSIRO-CCAM-r3355-r240x120-ISIMIP2b-AWAP')

#Main path
path <- args[1] #"/g/data/w35/amu561/Steven_CABLE_runs"

#AWRA inputs
awra_path <- args[2] #"/g/data/wj02/COMPLIANT/HMINPUT/output/AUS-5/BoM/"

#output on scratch
out_path <- args[3] 


#Run to process
model      <- args[4]
experiment <- args[5]
bc_method  <- args[6]


#Source functions
source(paste0(path, "/scripts/R/functions/parallel_process_flt.R"))
source(paste0(path, "/scripts/R/functions/write_daily_flt.R"))
source(paste0(path, "/scripts/R/functions/Calculate_VPD.R"))
#source(paste0(path, "/scripts/functions/Calculate_LWdown.R"))



#Use these options to read AWRA inputs
# DOMAIN = AUS-5
# ENSMEMBER = r1i1p1
# VERSION = latest


#Extra variables that need to be calculated
#(weather generator takes VPD (hPa) as input, converts it to Qair)
#Need 9am and 3pm values
#(also need lwdown but calculate this from weather-generated values)
#extra_vars <- c("vpd_09", "vpd_15")



#####################
### Process files ###
#####################


#Initialise cores
cl <- makeCluster(getOption('cl.cores', 6))



### list models ###

models <- list.files(awra_path)

 
#Parallel export
clusterExport(cl, list('model', 'experiment', 'bc_method',
                       'parallel_process_flt',
                       'write_daily_flt_lwdown',
                       'write_daily_flt_vpd',
                       'write_daily_flt',
                       'calculate_vpd',
                       #'calculate_lwdown',
                       'awra_path'))

clusterEvalQ(cl, library(raster))



# #Loop through models
# for (m in 1:length(models)) {
#   
#   
#   ### Get experiments ###
#   
#   experiments <- list.files(paste0(awra_path, "/", models[m]))
# 
#   clusterExport(cl, 'experiments')
#   
#   
  # #Loop through experiments
  # for (e in 1:length(experiments)) {
  #   
  #   
  #   #### Get bias correction methods ###
  #   
  #   bc_methods <- list.files(paste(awra_path, models[m], experiments[e], "r1i1p1", sep="/"))
  #     
  #   clusterExport(cl, 'bc_methods')
  #   
    
    # #Loop through BC methods
    # for (b in 1:length(bc_methods)) {
    #   
      
#Output directory for all variables
outdir_all <- paste(out_path, "/CABLE_inputs/Weather_generator_inputs", model,
                    experiment, bc_method, sep="/")

#If folder exists, remove to avoid potential issues
if (dir.exists(outdir_all)) unlink(outdir_all, recursive=TRUE)

#Create directory
dir.create(outdir_all, recursive=TRUE)


#Progress
print(paste0("Processing model: ", model, ", experiment: ", experiment, 
             ", bc_method: ", bc_method))


### Get variables ###

variables <- list.files(paste(awra_path, model, experiment, "r1i1p1", 
                        bc_method, "latest/day", sep="/"))

#Add extra variables
#variables <- append(variables, extra_vars)


clusterExport(cl, list('variables', 'outdir_all'))


#Parallel process variables
parLapply(cl, variables, function(x) parallel_process_flt(variable=x, awra_path, 
                                                          model=model, 
                                                          experiment=experiment, 
                                                          bc_method=bc_method,
                                                          outdir_all=outdir_all))

      
      # 
      # parallel_process_flt(variable="vpd_09", awra_path,
      #                      model=models[m],
      #                      experiment=experiments[e],
      #                      bc_method=bc_methods[b],
      #                      outdir_all=outdir_all)

#          
#     } #bc methods
# 
#   } #experiments
#   
# } #models
# 


stopCluster(cl)



library(raster)
library(parallel)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/"


#Source functions
source(paste0(path,"/scripts/R/functions/areal_mean.R"))


#Set percentile and scale
percentile <- "Perc_15"
scale      <- 3

#Variables
vars <- c("pr", "qtot", "sm_root")

#List metrics
metrics <- c("duration", "rel_intensity")

#Experiments
exp <- c("rcp45", "rcp85")

#Output directory
outdir <- paste0(path, "/NRM_means/")
dir.create(outdir)

outdir_obs <- paste0(path, "/obs/")
dir.create(outdir_obs)

outdir_mod <- paste0(path, "/models/")
dir.create(outdir_mod)


###################
### NRM regions ###
###################

nrm_regions <- raster("/g/data/wj02/MISC/NRM_CLUSTERS/NRM_clusters.nc")

#Set 0 values to NA (used for ocean masking)
nrm_regions[nrm_regions == 0] <- NA

#Mask NRM regions (some ocean pixels included in some NRM regions)
#Mask obs data (use an AWRA input file for masking for simplicity)
mask_data <- raster(paste0(path, "/AWRA_fields/",
                           "AWRA_fractional_vegetation_cover_monthly_climatology_1960_2005.nc"))

nrm_regions <- mask(nrm_regions, mask_data)

#Get region values for extracting later
nrm_vals <- sort(unique(values(nrm_regions)))
nrm_vals <- nrm_vals[which(!is.na(nrm_vals))]

nrm_labels <- c("CentralSlopes", "EastCoast", 
                "MurrayBasin", "MonsoonalNorth",
                "Rangelands", "SouthernSlopes",
                "SouthernAndSouthwesternFlatlands",
                "WetTropics")





#Loop through metrics
for (m in 1:length(metrics)) {
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    ####################
    ### Observations ###
    ####################
    
    #Either AGCD observed rainfall for pr or AWRA reference run for qtot/sm
    if (vars[v] == "pr") {
      
      obs_data <- brick(paste0(path, "/drought_metrics_AGCD/", scale, "-month/",
                               "drought_metrics_AGCD_precip_1900_2021_baseline_1970_2005_scale_",
                               scale, ".nc"), varname=metrics[m])
      
      #Crop to start in 1970 and end in 2020 (to match AWRA reference run years)
      #netcdf years are wrong so setting indices manually (data runs 1900-2021)
      obs_data <- obs_data[[841:1452]]
      
    } else {
      
      obs_data <- brick(list.files(paste0(path, "/drought_metrics_AWRA_ref/", scale, "-month/"),
                                   pattern=paste0("/drought_metrics_AWRA_ref_", vars[v]), full.names=TRUE))
      
    }
 
    
    ### Calculate areal means for each region ###
    
    #Make sure obs_data and nrm_regions have the same extent, otherwise the below
    #code will fail
    if(!extent(obs_data) == extent(nrm_regions)) {
      #pr obs data has a larger extent, need to crop the obs rather than NRM
      #(even though it's slower)
      obs_data <- crop(obs_data, nrm_regions)
    }
    
    #Loop through regions
    for (r in 1:length(nrm_vals)) {
      
      region_temp_file_obs <- paste0(outdir_obs, "/", nrm_labels[r], "_", 
                                     vars[v], "_obs-ref_", metrics[m], 
                                     "_scale_", scale, "_obs.csv")
    
      if (!file.exists(region_temp_file_obs)) {
        
        #Calculate areal mean of each time step
        obs_region_data <- areal_mean(mask(obs_data, nrm_regions,
                                           maskvalue=r, updatevalue=NA,
                                           inverse=TRUE))
        #Save data
        write.csv(obs_region_data, region_temp_file_obs)
        
      }
      
    }
    
    
    ##############
    ### Models ###
    ##############
    
    #Loop through experiments
    for (e in 1:length(exp)) {
      
      print(paste0("e:", e))
      
      
      #################
      ### Load data ###
      #################
      
      
      #List BC methods
      bc_methods <- list.files(paste0(path, "/drought_metrics/", scale, "-month"))
      
      #Find files
      bc_files <- lapply(bc_methods, function(x) list.files(paste0(path, "/drought_metrics/", 
                                                                   scale, "-month/", x),
                                                            pattern=glob2rx(paste0("drought_metrics_*", 
                                                                                   vars[v], "_", exp[e], "_", scale, ".nc")),
                                                            full.names=TRUE, recursive=TRUE))
      
      #Load data
      bc_data <- lapply(bc_files, function(x) lapply(x, function(f) brick(f, varname=metrics[m])))
      
      # 
      # #Crop bc data to run from 1970 to 2020 to match obs/reference data
      # bc_data <- lapply(bc_data, function(x) lapply(x, function(y) y[[121:732]]))
      # 
      # if(nlayers(obs_data) != nlayers(bc_data[[1]][[1]])) {
      #   stop("model and obs layers don't match")
      # }
      
      
      ### Loop through regions ###
      
      for (r in 1:length(nrm_vals)) {
        
        
        ### Model data ###
        
        #calculate areal means
        #Save output to speed up code
        region_temp_file <- paste0(outdir_mod, "/", nrm_labels[r], "_", 
                                   vars[v], "_", exp[e], "_", metrics[m], 
                                   "_scale_", scale, ".rds")
        
        if (!file.exists(region_temp_file)) {
          
          print(paste0("r:", r))

          #Parallelise
          cl <- makeCluster(getOption('cl.cores', 4))

          clusterExport(cl, c('nrm_regions', 'r', 'areal_mean', 'nlayers'))
          clusterEvalQ(cl, library(raster))

      
            # region_data <- lapply(bc_data, function(bc) lapply(bc, function(gcm) sapply(1:nlyr(gcm), 
          #                                                       function(n) areal_mean(mask(gcm[[n]], nrm_regions,
          #                                                                                         maskvalue=r,
          #                                                                                         updatevalue=NA,
          #                                                                                         inverse=TRUE)))))
          
          
          
          #Calculate areal meansmask(test, nrm_regions, maskvalue=NA, updatevalue=NA, inverse=TRUE)
          # region_data <- parLapply(cl=cl, bc_data, function(bc) lapply(bc, function(gcm) areal_mean(mask(gcm, nrm_regions,
          #                                                                                                maskvalue=r, updatevalue=NA,
          #                                                                                                inverse=TRUE))))
          #Calculate areal meansmask(test, nrm_regions, maskvalue=NA, updatevalue=NA, inverse=TRUE)
          region_data <- parLapply(cl=cl, bc_data, function(bc) lapply(bc, function(gcm) areal_mean(mask(gcm, nrm_regions,
                                                                                               maskvalue=r, updatevalue=NA,
                                                                                               inverse=TRUE))))
          
          #           # This works, but cannot be run from Rstudio
          # plan(multicore, workers=2)
          # future_lapply(s_list, FUN = function(x) sapply(1:nlayers(x), terra::global(data[[x]], mean, weights=cellSizes(data[[x]]))))
          # plan(sequential) # unsets the future cluster
          # 
          # 
          
          names(region_data) <- bc_methods
          
          #Save data
          saveRDS(region_data, region_temp_file)
          
          stopCluster(cl)
          
        } 
        
      } #regions
    } #experiments
  } #variables
} #metrics


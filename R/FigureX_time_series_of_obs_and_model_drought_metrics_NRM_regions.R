library(raster)
library(RColorBrewer)
library(maptools)
library(maps)
library(zoo)
library(parallel)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"


#Source functions
source(paste0(path,"/scripts/R/functions/areal_mean.R"))


#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3


#Variables
vars <- c("pr")#, "qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")

unit <- c("months", "% points") # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))


#Experiments
exp <- c("rcp45", "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Output directory
outdir <- paste0(path, "/Figures")
dir.create(outdir)




###################
### NRM regions ###
###################

nrm_regions <- raster("/g/data/wj02/MISC/NRM_CLUSTERS/NRM_clusters.nc")

#Set 0 values to NA (used for ocean masking)
nrm_regions[nrm_regions == 0] <- NA

#Get region values for extracting later
nrm_vals <- sort(unique(values(nrm_regions)))
nrm_vals <- nrm_vals[which(!is.na(nrm_vals))]

nrm_labels <- c("Central Slopes", "East Coast", 
                "Murray Basin", "Monsoonal North",
                "Rangelands", "Southern Slopes",
                "Southern and South-Western Flatlands",
                "Wet Tropics")

#cols=c("red", "blue", "purple", "green", "brown", "yellow", "orange", "black")
#plot(nrm_regions, breaks=c(0.5,1.2, 2.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), col=cols)


#Plotting colours
bc_cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")





#Save temporary outputs to speed up code
outdir_temp_obs <- paste0(path, "/temp/")

outdir_temp <- paste0(path, "/temp_NRM/")
dir.create(outdir_temp)


#Loop through metrics
for (m in 1:length(metrics)) {
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    #########################
    ### Load observations ###
    #########################
    
    
    obs_data <- brick(paste0(path, "/drought_metrics_AGCD/", scale, "-month/",
                            "drought_metrics_AGCD_precip_1900_2021_baseline_1970_2005_scale_",
                            scale, ".nc"), varname=metrics[m])
    
    
    
    #Crop to start in 1970
    #netcdf years are wrong so setting indices manually (data runs 1900-2021)
    obs_data <- obs_data[[841:nlayers(obs_data)]]
    
    #Mask obs data (use an AWRA input file for masking for simplicity)
    mask_data <- raster(paste0(path, "/AWRA_fields/",
                               "AWRA_fractional_vegetation_cover_monthly_climatology_1960_2005.nc"))
    
    
    obs_mask_file <- paste0(outdir_temp_obs, "/obs_masked_", vars[v], ".nc")
    
    if (!file.exists(obs_mask_file)) {
      #(need to first crop AGCD data as slightly larger extent than AWAP runs)
      obs_data_masked <- mask(crop(obs_data, mask_data), mask_data)
      
      writeRaster(obs_data_masked, obs_mask_file, overwrite=TRUE)

    } else {
      obs_data_masked <- brick(obs_mask_file)
    }
    
   
    ### Set up figure ###
    png(paste0(outdir, "/FigureX_NRM_time_series_of_historical_", metrics[m], 
               "_models_and_obs_perc_", percentile, "_scale_", scale, "_", 
               vars[v], ".png"),
        height=8.5, width=8.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.5, 0.4, 0.1))
    
    par(mfcol=c(length(nrm_vals), length(exp)))
    
    
    
    #Loop through experiments
    for (e in 1:length(exp)) {
    
    
      # layout(matrix(c(1:4, 5, 6:9, 5, 10:13, 5), nrow=5), heights=c(1, 1, 1, 1, 0.3,
      #                                                               1, 1, 1, 1,
      #                                                               1, 1, 1, 1))
      # #par(mfcol=c(3, 3))
      # par(bty="n")
      
    
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
      
      
      #Mask precip with NRM regions (precip includes some ocean grid cells that don't get
      #masked out if don't do this)
      
      if (vars[v] == "pr") {
        
        pr_temp_file <- paste0(outdir_temp, "/Masked_temp_", vars[v], "_", 
                               exp[e], "_", metrics[m], "_scale_", scale, ".rds")
        
        
        if (!file.exists(pr_temp_file)) {
          
        #Parallelise
        cl <- makeCluster(getOption('cl.cores', 4))
        
        clusterExport(cl, c('bc_data', 'nrm_regions'))
        clusterEvalQ(cl, library(raster))  
      
        #Mask 
        bc_data <- parLapply(cl=cl, bc_data, function(x) lapply(x, function(y) mask(y, nrm_regions)))
        
        # 
        # ### Raster ###
        # library(raster)
        # bc_data <- lapply(bc_files, function(x) lapply(x, function(f) brick(f, varname=metrics[m])))
        # #Parallelise
        # cl <- makeCluster(getOption('cl.cores', 4))
        # 
        # clusterExport(cl, c('nrm_regions'))
        # clusterEvalQ(cl, library(raster))  
        # 
        # #Mask 
        # bc_data <- parLapply(cl=cl, bc_data, function(x) lapply(x, function(y) raster::mask(y, nrm_regions)))
        # 
        # 
        # ### Terra ###
        # library(terra)
        # bc_data <- lapply(bc_files, function(x) lapply(x, function(f) rast(f, subds=metrics[m])))
        # #Parallelise
        # cl <- makeCluster(getOption('cl.cores', 4))
        # 
        # clusterExport(cl, c('nrm_regions'))
        # clusterEvalQ(cl, library(terra))  
        # 
        # #Mask 
        # bc_data <- parLapply(cl=cl, bc_data, function(x) lapply(x, function(y) terra::mask(y, nrm_regions)))
        # 
        # 
        # 

        #Save
        saveRDS(bc_data, pr_temp_file)
        
        stopCluster(cl)
        
        } else {
          bc_data <- readRDS(pr_temp_file)
        }
        
        
      }
      
      
      
      ### Loop through regions ###
      
      for (r in 1:length(nrm_vals)) {
        
        
        
          
        ### Model data ###
        
        #calculate areal means
        #Save output to speed up code
        region_temp_file <- paste0(outdir_temp, "/", nrm_labels[r], "_", 
                                   vars[v], "_", exp[e], "_", metrics[m], 
                                   "_scale_", scale, ".rds")
        
        if (!file.exists(region_temp_file)) {
          
          

          #Parallelise
          cl <- makeCluster(getOption('cl.cores', 4))

          clusterExport(cl, c('nrm_regions', 'r', 'areal_mean', 'nlayers'))
          clusterEvalQ(cl, library(raster))


          # #Calculate mean for NRM regions (parallelised where loops through years of each dataset)
          #For some reason this parallelisation doesn't work and ca'nt be bothered to work it out now
          # region_data <- parLapply(cl=cl, bc_data, function(bc) lapply(bc, function(gcm) sapply(1:nlyr(gcm), 
          #                                                                                  function(n) areal_mean_terra(mask(gcm[[n]], nrm_regions,
          #                                                                                                               maskvalue=r,
          #                                                                                                               updatevalue=NA,
          #                                                                                                               inverse=TRUE)))))
          #   
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
          region_data <- lapply(bc_data, function(bc) lapply(bc, function(gcm) areal_mean(mask(gcm, nrm_regions,
                                                                                               maskvalue=r, updatevalue=NA,
                                                                                               inverse=TRUE))))

          #           # This works, but cannot be run from Rstudio
          # plan(multicore, workers=2)
          # future_lapply(s_list, FUN = function(x) sapply(1:nlayers(x), terra::global(data[[x]], mean, weights=cellSizes(data[[x]]))))
          # plan(sequential) # unsets the future cluster
          # 
          # 
          
          #Save data
          saveRDS(region_data, region_temp_file)
          
          #stopCluster(cl)
          
        } else {
          region_data <- readRDS(region_temp_file)
        }
                                            
          
        
        ### Obs data ###
        
        region_temp_file_obs <- paste0(outdir_temp, "/", nrm_labels[r], "_", 
                                   vars[v], "_", exp[e], "_", metrics[m], 
                                   "_scale_", scale, "_obs.rds")
        
        if (!file.exists(region_temp_file_obs)) {
    
          
          
          #test=global(obs_data_masked, "mean")
          
          obs_region_data <- areal_mean(mask(obs_data_masked, nrm_regions,
                                             maskvalue=r, updatevalue=NA,
                                             inverse=TRUE))

          
          # obs_region_data <- sapply(1:nlyr(obs_data_masked), function(n)
          #                           areal_mean_terra(mask(obs_data_masked[[n]], nrm_regions,
          #                           maskvalue=r, updatevalue=NA, inverse=TRUE)))
          # 
          
          saveRDS(obs_region_data, region_temp_file_obs)
          
          
        } else {
          obs_region_data <- readRDS(region_temp_file_obs)
        }
        
        
        
        
        
        
        
      #   ################
      #   ### Plotting ###
      #   ################
      #   
      #   #Rolling mean factor
      #   k_fac <- 12*10
      #   
      # 
      #   yrs=1970:2099
      # 
      #   #Calculate obs rollmean
      #   #Add trend line
      #   x_obs  <- 1:length(obs_region_data)
      #   lm_obs <- lm(obs_region_data ~ x_obs)
      # 
      #   
      #   #Initialise
      #   mean_data <- list()
      #   max_data  <- list()
      #   min_data  <- list()
      #   
      #   #Loop through bc methods
      #   for (k in 1:length(region_data)) {
      # 
      #     mean_data[[k]] <- vector()
      #     max_data[[k]] <- vector()
      #     min_data[[k]] <- vector()
      # 
      #     #Loop through each time step to calculate ensemble mean and range. Not sure how to do this in lapply/sapply
      #     for (i in 1:length(region_data[[k]][[1]])) {
      # 
      #       mean_data[[k]][i] <- mean(sapply(region_data[[k]], function(x) x[i]), na.rm=TRUE)
      #       max_data[[k]][i]  <- max(sapply(region_data[[k]], function(x) x[i]), na.rm=TRUE)
      #       min_data[[k]][i]  <- min(sapply(region_data[[k]], function(x) x[i]), na.rm=TRUE)
      # 
      #     }
      #   }
      #  
      #   
      #   #Work out plot range
      #   plot_range <- range(c(rollmean(obs_region_data, k=k_fac),
      #                         unlist(sapply(mean_data, rollmean, k=k_fac))),
      #                       na.rm=TRUE)
      #   
      #   #Plot observations
      #   
      #   plot(rollmean(obs_region_data, k=k_fac, fill=NA, na.rm=TRUE), type='l',
      #        xlim=c(1, length(yrs)*12), ylim=plot_range,
      #        col="black", lty=1, lwd=2,
      #        ylab=paste0(metrics[m], " (", unit[m], ")"))
      #   
      #   #add NA as predict() ignores NA values produced in mean-vector
      #   lines(lm_obs$model$x_obs,  predict(lm_obs), col="black", lty=2, lwd=2 )
      #   
      #   
      #   #Loop through bc methods to plot
      #   for (k in 1:length(region_data)) {
      #     
      #     #plot mean and range
      #     lines(rollmean(mean_data[[k]], k=k_fac, fill=NA, na.rm=TRUE), col=bc_cols[k])
      # 
      #     #Calculate and plot trend
      #     x <- 1:length(mean_data[[k]])
      #     lm_model <- lm(mean_data[[k]] ~ x)
      # 
      #     #add NA as predict() ignores NA values produced in mean-vector
      #     lines(lm_model$model$x,  predict(lm_model), col=bc_cols[k], lty=2, lwd=2)
      # 
      #   }
      # 
      #   
      #   #Add area label
      #   if (e==1) mtext(side=2, nrm_labels[r], line=3)
      #   
      #   #Experiment label
      #   if(r==1) mtext(side=3, exp[e], line=2)
      # 
      } #regions
      
    } #experiments
    
    dev.off()
  } #variables
} #metrics


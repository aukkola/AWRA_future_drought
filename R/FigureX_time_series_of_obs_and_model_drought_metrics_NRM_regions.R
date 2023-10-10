library(ncdf4)
library(RColorBrewer)
library(maptools)
library(maps)
library(zoo)
library(raster)
library(grDevices)


#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"


#Source functions
source(paste0(path,"/scripts/R/functions/load_ncdf_var.R"))
source(paste0(path,"/scripts/R/functions/PolygonNA.R"))


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
bc_cols <- c("#33a02c",  "yellow", "#1f78b4", "#6a3d9a")




#Loop through metrics
for (m in 1:length(metrics)) {
  
  #Loop through variables
  for (v in 1:length(vars)) {
    

    ### Set up figure ###
    png(paste0(outdir, "/FigureX_NRM_time_series_of_historical_", metrics[m], 
               "_models_and_obs_perc_", percentile, "_scale_", scale, "_", 
               vars[v], ".png"),
        height=9.5, width=8.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.5, 0.4, 0.1))
    
    par(mfcol=c(length(nrm_vals), length(exp)))
    
    
    
    #Loop through experiments
    for (e in 1:length(exp)) {
    
      
      # bc_methods <- list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
      #                                 metrics[m], "/", exp[e], "/", vars[v]))
      
      #GCMs
      gcms <- list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
                                metrics[m], "/", exp[e], "/", vars[v], "/CCAM/"))
      
      
      
      #Loop through regions
      for (r in nrm_vals) {
        
        ### Model data ###
         
        # bc_files <- lapply(bc_methods, function(x) list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
        #                                                              metrics[m], "/", exp[e], "/", vars[v], "/", x),
        #                                                       recursive=TRUE, full.names=TRUE, pattern=paste0("region", r, ".nc")))
        
        #The "+." in the pattern allows two patterns to be specified (in this case GCM and region number)       
        gcm_files <- lapply(gcms, function(x) list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
                                                                      metrics[m], "/", exp[e], "/", vars[v]),
                                                               recursive=TRUE, full.names=TRUE, pattern=paste0(x, ".+", "region", r, ".nc")))

        
        #Check that found the right number of files (should be 4, one for each GCM)
        #if(any(sapply(bc_files, length) != 4)) stop("wrong number of files")
        if(any(sapply(gcm_files, length) != 4)) stop("wrong number of model files")
        
        
        model_data <- lapply(gcm_files, function(x) as.data.frame(lapply(x, function(y) load_ncdf_var(y, metrics[m]))))
        
        
        ### Obs data ###
        
        #Use AWAP for obs
        if (vars[v] == "pr") {
          
          obs_file <- list.files(paste0(path, "/NRM_means/data_obs/scale_", scale, 
                                        "/", metrics[m], "/", vars[v]),
                                        pattern=paste0("region", r, ".nc"), full.names=TRUE)
          
        #Use AWRA reference runs for runoff/SM
        } else {
          
          
        }
        
        #Sanity check
        if (length(obs_file) != 1) stop("wrong number of obs files")
        
        obs_data <- load_ncdf_var(obs_file, metrics[m])
        
        #Crop precip to start in 1970 to match models (starts 1900)
        if (vars[v] == "pr") {
          
          obs_data <- obs_data[841:length(obs_data)]
        }
        
        
        
        
        ################
        ### Plotting ###
        ################
      
        #Rolling mean factor
        k_fac <- 12*10
      
        yrs=1970:2099
      
        
        
        ### Model data ###
         
        #Calculate ensemble mean and range
        mean_data <- lapply(model_data, rowMeans) #only calculate mean if all models available
        max_data  <- lapply(model_data, function(x) apply(x, MARGIN=1, max, na.rm=TRUE))
        min_data  <- lapply(model_data, function(x) apply(x, MARGIN=1, min, na.rm=TRUE))
         
        #Replace Inf values with NA and linearly interpolate NA values (otherwise plotting is hard)
        mean_data <- lapply(mean_data, function(x) na.approx(replace(x, which(is.infinite(x)), NA), na.rm=FALSE)) 
        max_data  <- lapply(max_data, function(x) na.approx(replace(x, which(is.infinite(x)), NA), na.rm=FALSE)) 
        min_data  <- lapply(min_data, function(x) na.approx(replace(x, which(is.infinite(x)), NA), na.rm=FALSE))
                                                                     
        #Calculate rolling means
        mean_data <- lapply(mean_data, function(x) rollapply(x, k_fac, FUN=mean, na.rm=TRUE))
        max_data  <- lapply(max_data, function(x) rollapply(x, k_fac, FUN=mean, na.rm=TRUE))
        min_data  <- lapply(min_data, function(x) rollapply(x, k_fac, FUN=mean, na.rm=TRUE))
         
         
        ### Obs data ###
         
        obs_region_data <- rollapply(obs_data, k_fac, FUN=mean, na.rm=TRUE)
         
           
        #Work out plot range
        plot_y_range <- range(c(obs_region_data, unlist(max_data), unlist(min_data)), #unlist(mean_data)), #unlist(max_data), unlist(min_data)),
                               na.rm=TRUE)
      
        x_mod <- 1:length(mean_data[[1]])
        plot_x_range <- range(x_mod)
        
  
        ### Plot observations ###
      
        #Add trend line
        x_obs  <- 1:length(obs_region_data)
        lm_obs <- lm(obs_region_data ~ x_obs)
        
        plot(x_obs, obs_region_data, type='l',
             xlim=plot_x_range, ylim=plot_y_range, xaxt="n",
             col="black", lty=1, lwd=2,
             ylab=paste0(metrics[m], " (", unit[m], ")"))
      
         
      
      
        ### Plot model data ### 
        
        
        #Loop through GCMs to plot
        for (k in 1:length(mean_data)) {
      
          # #Plot ranges
          PolygonNA(x_mod, upper=max_data[[k]], lower=min_data[[k]], col=adjustcolor(bc_cols[k], alpha.f=0.3))

# 
#           if(all(is.na(mean_data[[k]]))) {
#             print(paste0("exp: ", e, ", model: ", k, ", region: ", r, ", all NA values"))
#             
#           } else {
#             
#             #plot mean and range
#             lines(mean_data[[k]], col=bc_cols[k])
#             
#             
#             #lines(rollmean(min_data[[k]], k=k_fac, fill=NA, na.rm=TRUE), col="red")
#             #lines(rollmean(max_data[[k]], k=k_fac, fill=NA, na.rm=TRUE), col="green")
#             
#             
  

#           }
#       
        }
        
        
      #Then loop over models again to aadd trend lines (do this afterwards so they appear on top)
      for (k in 1:length(mean_data)) {
        
        #Calculate and plot trend
        lm_model <- lm(mean_data[[k]] ~ x_mod)
        
        #add NA as predict() ignores NA values produced in mean-vector
        lines(lm_model$model$x,  predict(lm_model), col=bc_cols[k], lty=1, lwd=2)
        
      } 
      
      #Plot obs again so on top
      lines(x_obs, obs_region_data, type='l', col="black", lwd=2)
      
      #add NA as predict() ignores NA values produced in mean-vector
      lines(lm_obs$model$x_obs,  predict(lm_obs), col="black", lty=2, lwd=2 )
        
        
       #Add area label
       if (e==1) mtext(side=2, nrm_labels[r], line=3)
    
       #Experiment label
       if(r==1) mtext(side=3, exp[e], line=2)
    
      } #regions
      
    } #experiments
    
    dev.off()
  } #variables
} #metrics


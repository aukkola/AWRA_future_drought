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
exp <- c("rcp45") #, "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Output directory
outdir <- paste0(path, "/Figures")
dir.create(outdir)


#Tidier GCM labels (not included the full name
#because "0" and "-" cause issues. Dealing with this later)
gcm_labels <- c(CNRM  = "CNRM-CM5",
                CSIRO = "ACCESS1-0",
                MIROC = "MIROC5",         
                NOAA  = "GFDL-ESM2M")

#Fix y-axis across all plots to allow comparison of magnitudes
y_ranges_hist  <- list(duration = c(1.7, 2.6),
                       rel_intensity = c(40, 80))
y_ranges_trend <- list(duration = c(-0.001, 0.001),
                       rel_intensity = c(-0.006, 0.006))
  

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

#Symbols for model points
symbols <- c(0, 1, 2, 5)


#Loop through metrics
for (m in 1:length(metrics)) {
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    ### Set up figure ###
    png(paste0(outdir, "/FigureX_NRM_time_series_of_historical_", metrics[m], 
               "_models_and_obs_perc_", percentile, "_scale_", scale, "_", 
               vars[v], ".png"),
        height=11.7, width=8.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.5, 0.1))
    par(omi=c(0.1, 0.5, 0.2, 0.1))
    
    par(mfcol=c(4, 2)) #length(nrm_vals)/2, length(exp)))
    
    
    
    #Loop through experiments
    for (e in 1:length(exp)) {
      
      #Get names of bc methods for plotting
      bc_methods <- list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
                                       metrics[m], "/", exp[e], "/", vars[v]))
      
      #GCMs
      gcms <- list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
                                metrics[m], "/", exp[e], "/", vars[v], "/CCAM/"))
      
      
      
      #Loop through regions
      for (r in 1:length(nrm_vals)) {
        
        ### Model data ###
        
        # bc_files <- lapply(bc_methods, function(x) list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
        #                                                              metrics[m], "/", exp[e], "/", vars[v], "/", x),
        #                                                       recursive=TRUE, full.names=TRUE, pattern=paste0("region", r, ".nc")))
        
        #The "+." in the pattern allows two patterns to be specified (in this case GCM and region number)       
        gcm_files <- lapply(gcms, function(x) list.files(paste0(path, "/NRM_means/data/scale_", scale, "/",
                                                                metrics[m], "/", exp[e], "/", vars[v]),
                                                         recursive=TRUE, full.names=TRUE, pattern=paste0(x, ".+", "region", nrm_vals[r], ".nc")))
        
        
        #Check that found the right number of files (should be 4, one for each GCM)
        #if(any(sapply(bc_files, length) != 4)) stop("wrong number of files")
        if(any(sapply(gcm_files, length) != 4)) stop("wrong number of model files")
        
        
        model_data <- lapply(gcm_files, function(x) lapply(x, function(y) load_ncdf_var(y, metrics[m])))
        
        
        ### Obs data ###
        
        #Use AWAP for obs
        if (vars[v] == "pr") {
          
          obs_file <- list.files(paste0(path, "/NRM_means/data_obs/scale_", scale, 
                                        "/", metrics[m], "/", vars[v]),
                                 pattern=paste0("region", nrm_vals[r], ".nc"), full.names=TRUE)
          
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
 
        
        ### Model data ###
        x_mod <- 1:length(model_data[[1]][[1]])
        
        #Calculate trend over all years
        #Save trend, its error estimate and p-value
        model_trends <- lapply(model_data, function(x) lapply(x, function(y) summary(lm(y ~ x_mod))$coefficients[c(2,4,8)]))
          
        
        #Calculate historical mean (1970-2005)
        hist_indices <- 1:432
          
        model_hist_mean <- lapply(model_data, function(x) unlist(lapply(x, function(y) mean(y[hist_indices], na.rm=TRUE))))
        
        #Obs mean
        obs_mean <- mean(obs_data[hist_indices], na.rm=TRUE)
        
        
        
        #Initialise
        
        #hist_range <- range(c(obs_mean, unlist(model_hist_mean)))
        
        plot(c(0.5, 9.5), y_ranges_hist[[metrics[m]]], xaxt="n", 
             ylab="", xlab="", type="n")
        
        #Add historical mean line
        lines(c(-1, 5), rep(obs_mean, 2), lty=2)
        
        #Add vertical dividing line
        lines(c(5,5), c(100000, -100000))
        
        
        ### Plot historical ###
        
        #Model bars
        for(k in 1:length(model_hist_mean)) {
          
          polygon(c(k-0.35, k+0.35, k+0.35, k-0.35),
                  c(rep(max(model_hist_mean[[k]]), 2), rep(min(model_hist_mean[[k]]), 2)),
                    col=adjustcolor(bc_cols[k], alpha.f=0.7))
          
          #Add individual simulations
          points(rep(k, length(model_hist_mean[[k]])),  model_hist_mean[[k]],
                     pch=symbols, cex=1.5)
          
        }
        
        ### Plot trend ###
        
        ## Allow a second plot on the same graph
        par(new=TRUE)
        
        #Use x-axis from 5 to 9
        offset <- 5
        
        #trend_range <- range(unlist(lapply(model_trends, function(x) lapply(x, function(y) y[1]))))

        plot(c(0.5, 9.5), xaxt="n", yaxt="n", y_ranges_trend[[metrics[m]]],
             ylab="", xlab="", type="n", axes=FALSE)
        
        
        axis(side=4)
        
        #Add zero line
        lines(c(offset, 10), rep(0, 2), lty=2)
        
      
        #Model bars
        for(k in 1:length(model_hist_mean)) {
          
          gcm_vals <- unlist(lapply(model_trends[[k]], function(x) x[1]))
          pvals <- unlist(lapply(model_trends[[k]], function(x) x[3]))
          
          polygon(c(k-0.35+offset, k+0.35+offset, k+0.35+offset, k-0.35+offset),
                  c(rep(max(gcm_vals), 2), rep(min(gcm_vals), 2)),
                  col=bc_cols[k]) #adjustcolor(bc_cols[k], alpha.f=0.7))
          
          all_symbols <- symbols
          sign_symbols <- c(15, 16, 17, 23)
          
          if (any(pvals <= 0.05)) {
            all_symbols[which(pvals <= 0.05)] <- sign_symbols[which(pvals <= 0.05)]
          }
          
          #Add individual simulations
          points(rep(k+offset, length(gcm_vals)),  gcm_vals,
                 pch=all_symbols, bg="black", cex=1.5)
          
        }
        
        if(r==1) legend("topleft", legend=bc_methods, pch=symbols)
        
        
        #Region label
        mtext(side=3, nrm_labels[r], line=0, xpd=NA)
        
      } #regions
      
    } #experiments
    
    #Add x-axis (only bottom figure)
    axis(side=1, labels=rep(gcm_labels, 2), at=c(1:4, 6:9))
    
    
    dev.off()
  } #variables
} #metrics


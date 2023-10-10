library(raster)
library(RColorBrewer)
library(maptools)
library(maps)
library(parallel)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"

source(paste0(path, "/scripts/R/functions/pdf_and_skill.R"))


#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3


#Variables
vars <- c("pr", "qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")


#Experiments
exp <- c("rcp45")#, "rcp85")

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



#####################
### Plot settings ###
#####################


cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")


unit <- c("months", "% points") # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7

#Can only specify one experiment for this code
if(length(exp) > 1) stop("too many experiments")

#Tidier GCM labels
gcm_labels <- c(CNRM  = "CNRM-CM5",
                CSIRO = "ACCESS1-0",
                MIROC = "MIROC5",         
                NOAA  = "GFDL-ESM2M")


#Set up parallel processing for trend calculation
cl <- makeCluster(getOption('cl.cores', 16))

clusterEvalQ(cl, library(raster))


#Plot GCM and BC in separate plot
plot_type <- c("GCM", "BC")


#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  #Progress
  print(paste0("plotting metric ", m, "/", length(metrics)))

  #Plot type (GCM/BC)  
  for (p in 1:length(plot_type)) {
  

    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_historical_PDFs_of_", plot_type[p], "_in_", metrics[m], "_",
               percentile, "_", scale, "_", exp, "_by_NRM_region.png"),
        height=9.5, width=6.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.4, 0.4, 0.3))
    
    par(mfrow=c(length(vars), 2))
    
    # layout(matrix(c(1:4, 5, 6:9, 5, 10:13, 5), nrow=5), heights=c(1 ,0.3, 1, 1, 0.3,
    #                                                               1 ,0.3, 1, 1, 0.3,
    #                                                               1 ,0.3, 1, 1, 0.3))
    # #par(mfcol=c(3, 3))
    #par(bty="n")
    
    
    
    #Loop through variables
    for (v in 1:length(vars)) {
      
      ################
      ### Obs data ###
      ################
      
      if (vars[v] == "pr") {
        
        obs_file <- paste0(path, "/drought_metrics_AGCD/", scale, "-month/",
                           "drought_metrics_AGCD_precip_1900_2021_baseline_1970_2005_scale_3.nc")
        
        #Select years 1970-2020 to match model runs and AWRA reference runs
        obs_data <- brick(obs_file, varname=metrics[m])[[841:1452]]
        
      } else {
        
        obs_file <- paste0(path, "/drought_metrics_AWRA_ref/", scale, "-month/",
                           "drought_metrics_AWRA_ref_", vars[v], "_scale_", scale, "_1960_2020.nc")
        obs_data <- brick(obs_file, varname=metrics[m])
      }
      
      #Calculate obs mean
      obs_mean <- mean(obs_data, na.rm=TRUE)
      
      
      
      ##################
      ### Model data ###
      ##################
      
      data_files <- list.files(paste0(path, "/drought_metrics/", 
                                      scale, "-month/"),
                               pattern=paste0(vars[v], "_", exp),
                               full.names=TRUE, recursive=TRUE)
      
      #Should find 16 files in each case, check
      if (length(data_files) != 16) {
        stop("Incorrect number of files found")
      }
      
      #Get BC methods
      bc_methods <- list.files(paste0(path, "/drought_metrics/", 
                                      scale, "-month/"))
      
      #Get GCMs
      gcms <- list.files(paste0(path, "/drought_metrics/", 
                                scale, "-month/CCAM/"))
      
      
      #Load data
      #Only get data for 1970- to match observations
      data <- lapply(data_files, function(x) brick(x, varname=metrics[m]))
      
      
      #Calculate model means
      model_mean <- parLapply(cl, data, mean, na.rm=TRUE)
      
      
      #Mask obs with models (obs includes ocean areas)
      obs_mean <- mask(crop(obs_mean, model_mean[[1]]),  model_mean[[1]])
      
      
      
      ################
      ### Plotting ###
      ################
      
      #Get area of each pixel
      area <- values(mask(area(obs_mean), obs_mean))
      
      #Get data values
      obs_vals <- values(obs_mean)
      mod_vals <- lapply(model_mean, values)
      
      
      #Set bins for frequency calculation
      nbins <- 50
      
      data_range <- range(obs_vals, unlist(mod_vals), na.rm=TRUE)
      
      breaks <- seq(data_range[1], data_range[2], length.out=nbins)
      
      #Calculate obs density
      obs_freqs <- freqs(obs_vals, breaks, area)
      
      # 
      # #Calculate densities (fraction of land area in each bin)
      # mod_freqs <- lapply(mod_vals, function(x) freqs(x, breaks, area))
      # 
      # 
      # 
      # #Calculate Perkins skill score
      # perkins_skill <- lapply(mod_freqs, function(x) perkins_skill_score(x, obs_freqs)) 
      # 
      # #Calculate normalised mean error
      # #Not using this as doesn't give the sign of bias
      # #nme <- lapply(mod_vals, function(x) NME(x, obs_vals))
      # 
      # #Using mean bias error for now
      # mbe <- lapply(mod_vals, function(x) mean(x, na.rm=TRUE) - mean(obs_vals, na.rm=TRUE))
      # 
      # 
      
      #######################
      ### Plot PDF by GCM ###
      #######################
      
      gcm_ind <- lapply(gcms, function(x) which(grepl(x, data_files)))
      
      gcm_freqs <- lapply(gcm_ind, function(x) freqs(unlist(mod_vals[x]), breaks, rep(area, length(x))))
      
      
      #Get plot ranges
      plot_y_range <- range(unlist(gcm_freqs), obs_freqs, na.rm=TRUE)
      
      #x-range has some crazy values from the percentage calculation, ignore these
      #by taking x-range as the 10/90th percentile
      plot_x_range <- c(1, scale+3) #range(breaks) #quantile(unlist(lapply(hist_gcm, function(x) x$x)), 
      #probs=c(0.1, 0.9), na.rm=TRUE)
      
      
      #Plotting
      plot(plot_x_range, plot_y_range, type="n", ylab="n", xlab="n")
      
      #Work out x values (bin centres)
      offset <- breaks[2] - breaks[1]
      
      x_centred <- breaks[2:length(breaks)] - offset
      
      
      #Add model lines
      for (gcm in 1:length(gcm_freqs)) {
        
        #Plot line
        lines(x_centred, gcm_freqs[[gcm]], col=cols[gcm])
      }
      
      if(v==1) {
        mtext(side=3, line=1, cex=1.2, "GCM")
        
      }
      
      mtext(side=2, line=2, var_labels[v])
      
      
      #Add obs line
      lines(x_centred, obs_freqs, col="black", lwd=2)
      
      #Add legend with skill scores
      perkins_gcm <- sapply(gcm_freqs, function(x) perkins_skill_score(x, obs_freqs))
      mbe_gcm     <- sapply(gcm_ind, function(x) mean(unlist(mod_vals[x]), na.rm=TRUE) - mean(obs_vals, na.rm=TRUE))
      
      
      legend("topright", legend=paste0(gcm_labels, " (s=", round(perkins_gcm, digits=2), ", MBE=", 
                                       round(mbe_gcm, digits=2), ")"), col=cols, bty="n", lty=1, cex=0.8)
      
      
      #############################
      ### Plot PDF by BC method ###
      #############################
      
      bc_ind <- lapply(bc_methods, function(x) which(grepl(x, data_files)))
      
      #Calculate densities
      bc_freqs <- lapply(bc_ind, function(x) freqs(unlist(mod_vals[x]), breaks, rep(area, length(x))))
      
      #Get plot ranges
      plot_y_range <- range(unlist(lapply(bc_freqs, function(x) x$y)), na.rm=TRUE)
      
      #x-range has some crazy values from the percentage calculation, ignore these
      #by taking x-range as the 10/90th percentile
      #plot_x_range <- range(breaks) #quantile(unlist(lapply(hist_bc, function(x) x$x)), 
      #probs=c(0.1, 0.9), na.rm=TRUE)
      
      
      #Plotting
      plot(plot_x_range, plot_y_range, type="n", ylab="n", xlab="n")
      
      
      #Add model lines
      for (bc in 1:length(bc_freqs)) {
        
        #Plot line
        lines(x_centred, bc_freqs[[bc]], col=cols[bc])
      }
      
      if(v==1) {
        mtext(side=3, line=1, "BC method")
      }
      
      #Add obs line
      lines(x_centred, obs_freqs, col="black", lwd=2)
      
      
      #Add legend with skill scores
      perkins_bc <- sapply(bc_freqs, function(x) perkins_skill_score(x, obs_freqs))
      mbe_bc     <- sapply(bc_ind, function(x) mean(unlist(mod_vals[x]), na.rm=TRUE) - mean(obs_vals, na.rm=TRUE))
      
      
      legend("topright", legend=paste0(bc_methods, " (s=", round(perkins_bc, digits=2), ", MBE=", 
                                       round(mbe_bc, digits=2), ")"), col=cols, bty="n", lty=1, cex=0.8)
      
      
    } #variables
    
    dev.off ()
  } #GCM/BC
} #metrics



stopCluster(cl) #need to put this here, otherwise temporary .grd files craeted by parallel code become unavailable



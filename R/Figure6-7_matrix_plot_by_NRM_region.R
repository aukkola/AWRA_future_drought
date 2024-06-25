library(ncdf4)
library(plot.matrix)
library(raster)
library(ggplot2)
library(ggpubr)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"

#source function

source(paste0(path, "/scripts/R/functions/load_ncdf_var.R"))
source(paste0(path, "/scripts/R/functions/trend_per_pixel.R"))
source(paste0(path, "/scripts/R/functions/add_raster_legend.R"))

#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3


#Variables
vars <- c("pr", "qtot")#, "sm_root")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity", "timing")#, "frequency")


#Experiments
exp <- c("rcp45")#, "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Output directory
outdir <- paste0(path, "/Figures")
dir.create(outdir)


#####################
### Plot settings ###
#####################


hist_cols <- colorRampPalette(rev(c("#b2182b", "#d6604d", "#f4a582", "#fddbc7",
                  "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac")))

trend_cols <- colorRampPalette(rev(c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3",
                "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e")))
  
  
  
#drought metric units
unit <- c(duration="months", rel_intensity="%", timing="%")

#Neater labels
metric_labels <-c(duration="Duration", rel_intensity="Relative intensity", timing="Area under drought")


#Can only specify one experiment for this code
if(length(exp) > 1) stop("too many experiments")

#Tidier GCM labels
gcm_labels <- c(CNRM  = "CNRM-CM5",
                CSIRO = "ACCESS1-0",
                MIROC = "MIROC5",         
                NOAA  = "GFDL-ESM2M")

#Set breaks
minn <- -10000
maxx <- 10000

breaks_hist <- list(duration=c(minn, seq(-0.5, 0.5, by=0.1), maxx),
                    rel_intensity=c(minn, seq(-25, 25, by=5), maxx),
                    timing=c(minn, seq(-25, 25, by=5), maxx))

breaks_trend <- list(pr=list(duration=c(minn, seq(-0.05, 0.05, by=0.01), maxx), # seq(-0.5, 0.5, by=0.1), maxx)/1000,
                             rel_intensity=c(minn, seq(-25, 25, by=5), maxx)/1000,
                             timing=c(minn, seq(-25, 25, by=5), maxx)/1000),
                    qtot=list(duration=c(minn, seq(-0.16, 0.16, by=0.04), maxx), #seq(-1.5, 1.5, by=0.3), maxx)/1000,
                              rel_intensity=c(minn, seq(-25, 25, by=5), maxx)/1000,
                              timing=c(minn, seq(-25, 25, by=5), maxx)/1000),
                    sm=list(duration=c(minn, seq(-1, 1, by=0.1), maxx)/1000,
                            rel_intensity=c(minn, seq(-25, 25, by=5), maxx)/1000,
                            timing=c(minn, seq(-25, 25, by=5), maxx)/1000))


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
                "S/SW Flatlands",
                "Wet Tropics")




#Loop through variables
for (v in 1:length(vars)) {
  
  #Loop through metrics
  for (m in 1:length(metrics)) {

    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_matrix_NRM_regions_of_gmcs_and_bc_methods_in_", vars[v], "_", 
               metrics[m], "_", percentile, "_", scale, "_", exp, ".png"),
        height=9.5, width=8.3, units="in", res=400)
    
    
    par(mai=c(0.8, 0.2, 0.2, 0.2))
    par(omi=c(0.1, 1.2, 0.4, 0.3))
    
    layout(matrix(c(1,2, 3, 3, 4, 5, 6, 6), ncol=2, byrow=TRUE), heights=c(1,0.3,1,0.3))
    

    #Set up matrix to collate data (4 GCMs or 4 BC methods)
    hist_gcm_matrix  <- matrix(nrow=length(nrm_vals), ncol=4)
    hist_bc_matrix   <- matrix(nrow=length(nrm_vals), ncol=4)

    trend_gcm_matrix <- matrix(nrow=length(nrm_vals), ncol=4)
    trend_bc_matrix  <- matrix(nrow=length(nrm_vals), ncol=4)
    
    
    #Loop through regions
    for(r in 1:length(nrm_vals)) {
      
      ################
      ### Obs data ###
      ################
      
      #Get obs file
      obs_file <- list.files(paste0(path, "/NRM_means/data_obs/scale_", scale, "/", metrics[m],
                                    "/", vars[v], "/"), pattern=paste0("region", nrm_vals[r]), full.names=TRUE)
      
      #Select years 1970-2020 to match model runs and AWRA reference runs
      obs_data <- load_ncdf_var(obs_file, metrics[m])
      
      #Need to convert timing to fractional area (input data is the number of pixels under drought
      #in each NRM region, divide by the total number of pixels in that region)
      if (metrics[m] == "timing") {
        obs_data <- obs_data / length(which(values(nrm_regions) == nrm_vals[r])) *100
        if(any(obs_data > 100)) stop("timing gone wrong")
      }
      
      #Indices for years 1970-2020 for data running 1960 onwards (i.e. everything 
      #except obs precip)
      hist_ind <- 121:732
      
      #Extract 1970-2020 for obs
      if(vars[v] == "pr") {
        #Calculate obs mean
        obs_mean <- mean(obs_data[841:1452], na.rm=TRUE)

      } else {
      #Calculate obs mean
      obs_mean <- mean(obs_data[hist_ind], na.rm=TRUE)

      }

      
      
      
      ##################
      ### Model data ###
      ##################
      
      file_path <- paste0(path, "/NRM_means/data/scale_", 
                          scale, "/", metrics[m], "/", exp, "/", vars[v])
      
      data_files <- list.files(file_path, pattern=paste0("region", nrm_vals[r]),
                               full.names=TRUE, recursive=TRUE)
      
      #Should find 16 files in each case, check
      if (length(data_files) != 16) {
        stop("Incorrect number of files found")
      }
      
      #Get BC methods
      bc_methods <- list.files(paste0(file_path))
      
      #Get GCMs
      gcms <- list.files(paste0(file_path, "/CCAM/"))
      
      
      #Load data
      #Only get data for 1970- to match observations
      data <- lapply(data_files, function(x) load_ncdf_var(x, metrics[m]))
      
      #Convert to fractional area
      if (metrics[m] == "timing") {
        data <- lapply(data, function(x) x / length(which(values(nrm_regions) == nrm_vals[r])) * 100)
        if(any(unlist(data) > 100)) stop("timing gone wrong")
      }
      
      
      #Calculate model historical means
      model_mean <- lapply(data, function(x) mean(x[hist_ind], na.rm=TRUE))
    
        
      #Calculate model trend (using 1970-2099)
      model_trend <- lapply(data, function(x) trend_per_pixel(x[hist_ind[1]:length(x)])[1])
      
      
      ####################
      ### Group by GCM ###
      ####################
      
      #Find indices for each GCM
      gcm_ind <- lapply(gcms, function(x) which(grepl(x, data_files)))
      
      #Calculate the mean of GCM means, then take the difference from obs
      #Calculating the mean for each model run separarely first, I think this is the best way to do it.
      #Other option would be to lump all the GCM time series together and then take the mean. But because
      #some model runs will have more/less drought events than others, I think this could bias results
      gcm_bias <- lapply(gcm_ind, function(x) mean(unlist(model_mean[x]), na.rm=TRUE) - obs_mean)
      
      #Again need to take trends individually, then average. I think this is ok? can't really lump together
      #and then take trend
      gcm_trend <- lapply(gcm_ind, function(x) mean(unlist(model_trend[x])))
      
      
      #Add to plotting matrices
      hist_gcm_matrix[r,]  <- unlist(gcm_bias)
      trend_gcm_matrix[r,] <- unlist(gcm_trend)
      
      
      ###################
      ### Group by BC ###
      ###################
      
      #Find indices for each GCM
      bc_ind <- lapply(bc_methods, function(x) which(grepl(x, data_files)))
      
      #historical bias
      bc_bias <- lapply(bc_ind, function(x) mean(unlist(model_mean[x]), na.rm=TRUE) - obs_mean)
      
      #trend
      bc_trend <- lapply(bc_ind, function(x) mean(unlist(model_trend[x])))
      
      #Add to plotting matrices
      hist_bc_matrix[r,]  <- unlist(bc_bias)
      trend_bc_matrix[r,] <- unlist(bc_trend)
      
    }
    
    
    
    
    
    ################
    ### Plotting ###
    ################
    
    ### GCM and BC bias (top row, two panels) ###
    
    plot_data <- list(hist_gcm_matrix, hist_bc_matrix)
    
    plot_labs <- list(as.vector(gcm_labels), bc_methods)
    
    breaks <- breaks_hist[[metrics[m]]]
    
    main_title <- c("GCM", "DS-BC method")
    
    for (p in 1:length(plot_data)) {
      
      #matrix plot
      image(t(plot_data[[p]]), breaks=breaks, col=hist_cols(length(breaks)-1), xaxt="n",
            yaxt="n")
      
      #x-axis
      axis(side=1, labels=plot_labs[[p]], las=2, 
           at=seq(0,1, length.out=length(plot_labs[[p]])))
      
      #y-axis
      if(p==1) axis(side=2, labels=nrm_labels, las=2,
             at=seq(0,1, length.out=length(nrm_labels)))
        
      mtext(side=3, main_title[p], line=1, cex=1.2)
    
    }
    
    
    ### Legend ###
    
    plot(1, type="n", xaxt="n", yaxt="n", bty="n", ylab="", xlab="") #empty plot
    add_raster_legend2(hist_cols(length(breaks)-1), breaks[2:(length(breaks)-1)],
                       plot_loc=c(0.2,0.8,-12.95,-11), main_title=paste0("Bias (", unit[metrics[m]], ")"),
                       spt.cex=1.5, title_fac=-0.5)
    
    
    
    ### GCM and BC trend (bottom row, two panels) ###
    
    plot_data <- list(trend_gcm_matrix, trend_bc_matrix)
    
    plot_labs <- list(as.vector(gcm_labels), bc_methods)
    
    breaks <- breaks_trend[[vars[v]]][[metrics[m]]]
    
      
    for (p in 1:length(plot_data)) {
      
      #If duration, changes units to months per decade for tidier labels
      #(breaks already set accordingly)
      if (metrics[m] == "duration") plot_data[[p]] <- plot_data[[p]] * 120
        
      #matrix plot
      image(t(plot_data[[p]]), breaks=breaks, col=trend_cols(length(breaks)-1), xaxt="n",
            yaxt="n")
      
      #x-axis
      axis(side=1, labels=plot_labs[[p]], las=2, 
           at=seq(0,1, length.out=length(plot_labs[[p]])))
      
      #y-axis
      if(p==1) axis(side=2, labels=nrm_labels, las=2,
                    at=seq(0,1, length.out=length(nrm_labels)))
      
    }
    
    ### Legend ###
    
    #If duration, change units to months per year for tidier labels
    #ALREADY DONE ABOVE when setting limits
    #if (metrics[m] == "duration") breaks <- round(breaks *120, digits=2)
    
    
    plot(1, type="n", xaxt="n", yaxt="n", bty="n", ylab="", xlab="") #empty plot
    add_raster_legend2(trend_cols(length(breaks)-1), breaks[2:(length(breaks)-1)],
                       plot_loc=c(0.2,0.8,-12.95,-11), main_title=paste0("Trend (", unit[metrics[m]], "/decade)"),
                       spt.cex=1.5, xpd=NA, title_fac=-0.5)
    
    dev.off ()
    
    } #metrics
  
} #variables








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
vars <- c("pr", "qtot", "sm_root")

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")

metric_labels <- c("Duration (months)", "Relative intensity (%)")#, "frequency") #labels for plotting


#Experiments
exp <- c("rcp45")#, "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Output directory
outdir <- paste0(path, "/Figures")
dir.create(outdir)



#####################
### Plot settings ###
#####################


cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
#c("#e41a1c", "#377eb8", "#984ea3", "#ff7f00")
#("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")


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


x_range <- list(duration=list(pr=c(1,3.5) + scale-1,
                              qtot=c(1,5)+ scale-1,
                              sm_root=c(1,6)+ scale-1),
                rel_intensity=list(pr=c(10,90),
                                   qtot=c(10,90),
                                   sm_root=c(10,90)))

#Set up parallel processing for trend calculation
#cl <- makeCluster(getOption('cl.cores', 28))

#clusterEvalQ(cl, library(raster))
beginCluster(28)

#Plot GCM and BC in separate plot
#plot_type <- c("GCM", "BC")


#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  #Progress
  print(paste0("plotting metric ", m, "/", length(metrics)))
  
  # #Plot type (GCM/BC)  
  # for (p in 1:length(plot_type)) {
  #   
    
    ### Set up figure ###
    png(paste0(outdir, "/Figure5", "_historical_PDFs_in_", metrics[m], "_",
               percentile, "_", scale, "_", exp, "_whole_continent.png"),
        height=9.5, width=7.3, units="in", res=400)
    
    
    par(mai=c(0.3, 0.1, 0.2, 0.1))
    par(omi=c(0.5, 1.0, 0.4, 0.3))
    
    par(mfrow=c(length(vars), 2))
    
    # layout(matrix(c(1:4, 5, 6:9, 5, 10:13, 5), nrow=5), heights=c(1 ,0.3, 1, 1, 0.3,
    #                                                               1 ,0.3, 1, 1, 0.3,
    #                                                               1 ,0.3, 1, 1, 0.3))
    # #par(mfcol=c(3, 3))
    #par(bty="n")
    
    
    
    #Loop through variables
    for (v in 1:length(vars)) {
      
      print(paste0("variable ", v, "/", length(vars)))
      
      
      ################
      ### Obs data ###
      ################
      
      if (vars[v] == "pr") {
        
        obs_file <- paste0(path, "/drought_metrics_AGCD/", scale, "-month/",
                           "drought_metrics_AGCD_precip_1900_2021_baseline_1970_2005_scale_3.nc")
        
        #Select years 1970-2020 to match model runs and AWRA reference runs
        obs_data <- brick(obs_file, varname=metrics[m])[[841:1452]]
        
      } else {
        if(vars[v]=="sm_root"){
          var="sm"
        } else {
          var=vars[v]
        }
        #quick bug fix before re-run ref file
        obs_file <- paste0(path, "/drought_metrics_AWRA_ref/", scale, "-month/",
                           "drought_metrics_AWRA_ref_", var, "_scale_", scale, "_1960_2020.nc")
                           #vars[v], "_scale_", scale, "_1960_2020.nc")
        
        #Select 1970-2020
        obs_data <- brick(obs_file, varname=metrics[m])
        obs_data <- obs_data[[121:nlayers(obs_data)]]
      }
      
      #Calculate obs mean
      obs_mean <- clusterR(obs_data, mean, args=list(na.rm=TRUE))
      
      
       
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
      #Only get data for 1970-2020 to match observations (model data runs 1960-2099)
      data <- lapply(data_files, function(x) brick(x, varname=metrics[m])[[121:732]])
      
      start=Sys.time()
      #Calculate model means
      model_mean <- lapply(data, function(x) clusterR(x, mean, args=list(na.rm=TRUE))) #parLapply(cl, data, mean, na.rm=TRUE)
      end=Sys.time()
      print(paste("time: ", end-start))
      
      
      #Mask rainfall obs with models (obs includes ocean areas)
      if (vars[v]=="pr") {
        obs_mean <- mask(crop(obs_mean, model_mean[[1]]),  model_mean[[1]])
      } 
      
      
      #Need to add scale-1 to duration so numbers reflect the scale
      #e.g. a one-month drought for scale=3 should be 3 months
      if (metrics[m] == "duration") {
        obs_mean <- obs_mean + (scale-1)
        model_mean <- lapply(model_mean, function(x) x + (scale-1))
      }
      
      
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
      
      if(metrics[m] == "duration") {
        
        data_range <- c(1,6) #range(obs_vals, unlist(mod_vals), na.rm=TRUE)
      } else if (metrics[m] == "rel_intensity") {
        data_range <- c(0,100) #range(obs_vals, unlist(mod_vals), na.rm=TRUE)
      } else {
        data_range <- range(obs_vals, unlist(mod_vals), na.rm=TRUE)
      }
      
      data_range <- x_range[[metrics[m]]][[vars[v]]]
      
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
      
      
      #############################
      ### Calculate frequencies ###
      #############################
      
      #GCM
      gcm_ind <- lapply(gcms, function(x) which(grepl(x, data_files)))
      
      gcm_freqs <- lapply(gcm_ind, function(x) freqs(unlist(mod_vals[x]), breaks, rep(area, length(x))))
      
      #BC
      bc_ind <- lapply(bc_methods, function(x) which(grepl(x, data_files)))
      
      bc_freqs <- lapply(bc_ind, function(x) freqs(unlist(mod_vals[x]), breaks, rep(area, length(x))))
      
      
      #Get plot ranges
      plot_y_range <- range(unlist(gcm_freqs), bc_freqs, obs_freqs, na.rm=TRUE)
      plot_x_range <- x_range[[metrics[m]]][[vars[v]]] 
      
      
      #######################
      ### Plot PDF by GCM ###
      #######################
    
      
      #Plotting
      plot(plot_x_range, plot_y_range, type="n", ylab="n", xlab="n")
      
      #Work out x values (bin centres)
      offset <- breaks[2] - breaks[1]
      
      x_centred <- breaks[2:length(breaks)] - offset
      
      #Add obs line
      lines(x_centred, obs_freqs, col="black", lwd=1)
      
      #Add model lines
      for (gcm in 1:length(gcm_freqs)) {
        
        #Plot line
        lines(x_centred, gcm_freqs[[gcm]], col=cols[gcm])
      }
      
      if(v==1) {
        mtext(side=3, line=1, cex=1.2, "GCM")
        
      }
      
      #Variable label
      mtext(side=2, line=5, cex=1.2, var_labels[v])
      
      #y-label
      mtext(side=2, line=2.5, "Fraction of land area (-)")
      
      #x-label
      if(v==3) {
        mtext(side=1, line=2.5, metric_labels[m])
      }
      
      
      #Add legend with skill scores
      perkins_gcm <- sapply(gcm_freqs, function(x) perkins_skill_score(x, obs_freqs))
      mbe_gcm     <- sapply(gcm_ind, function(x) mean(unlist(mod_vals[x]), na.rm=TRUE) - mean(obs_vals, na.rm=TRUE))
      
      #Need to place legend in a different spot depending on variable/metric
      position <- "topright"
      
      if (metrics[m] == "rel_intensity" ){ #& vars[v] %in% c("pr", "qtot")) {
        position <- "topleft"
      } 
      
      legend(position, legend=paste0(gcm_labels, " (s=", round(perkins_gcm, digits=2), ", MBE=", 
                                       round(mbe_gcm, digits=2), ")"), col=cols, bty="n", lty=1, cex=0.8)
      
      
      #############################
      ### Plot PDF by BC method ###
      #############################
      
        
      #x-range has some crazy values from the percentage calculation, ignore these
      #by taking x-range as the 10/90th percentile
      #plot_x_range <- range(breaks) #quantile(unlist(lapply(hist_bc, function(x) x$x)), 
      #probs=c(0.1, 0.9), na.rm=TRUE)
      
      
      #Plotting
      plot(plot_x_range, plot_y_range, type="n", ylab="n", xlab="n", yaxt="n")
      
      #Add obs line
      lines(x_centred, obs_freqs, col="black", lwd=1)
      
      #Add model lines
      for (bc in 1:length(bc_freqs)) {
        
        #Plot line
        lines(x_centred, bc_freqs[[bc]], col=cols[bc])
      }
      
      #Main title
      if(v==1) {
        mtext(side=3, line=1, cex=1.2, "DS-BC method")
      }
      
      #x-label
      if(v==3) {
        mtext(side=1, line=2.5, metric_labels[m])
      }
      
      #Add legend with skill scores
      perkins_bc <- sapply(bc_freqs, function(x) perkins_skill_score(x, obs_freqs))
      mbe_bc     <- sapply(bc_ind, function(x) mean(unlist(mod_vals[x]), na.rm=TRUE) - mean(obs_vals, na.rm=TRUE))
    
      legend(position, legend=paste0(bc_methods, " (s=", round(perkins_bc, digits=2), ", MBE=", 
                                       round(mbe_bc, digits=2), ")"), col=cols, bty="n", lty=1, cex=0.8)
      
      
    } #variables
    
    dev.off()
  # } #GCM/BC
} #metrics

endCluster()

#stopCluster(cl) #need to put this here, otherwise temporary .grd files created by parallel code become unavailable
# 
# "#AE3C60","#F3C33C", "#267778", "#82B4BB"
# "#6EC3C1", "#335120", "#9DCC5F", "#0D5F8A"
# "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"


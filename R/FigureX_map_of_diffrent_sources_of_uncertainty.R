library(raster)
library(RColorBrewer)
library(maptools)
library(maps)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"


#Source functions
source(paste0(path,"/scripts/R/functions/add_raster_legend.R"))
source(paste0(path,"/scripts/R/functions/areal_mean.R"))



#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3

#baseline   <- "1970-2005_vs_2064-2099"


#Variables
vars <- c("pr", "qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")


#Experiments
exp <- c("rcp45", "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")


outdir <- paste0(path, "/Figures")
dir.create(outdir)


#####################
### Plot settings ###
#####################


#Set plot colours

#Historical mean
cols <- colorRampPalette(c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4",
  "#1d91c0", "#225ea8","#0c2c84"))
  
  #colorRampPalette(c("#edf8fb", "#bfd3e6", "#9ebcda", "#8c96c6",
  #"#8c6bb1", "#88419d","#6e016b"))
  
  #colorRampPalette(c("#ffffb2", "#feb24c",
   #                       "#fc4e2a", "#b10026"))
                                         

lims <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  
lims <- c(0, 20, 40, 60, 80, 100)
panel_cex=0.7



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  ### Set up figure ###
  png(paste0(outdir, "/FigureX", "_Sources_of_uncertainty_in_", metrics[m], "_",
             percentile, "_", scale, ".png"),
      height=6.5, width=8.3, units="in", res=400)
  
  
  par(mai=c(0, 0.1, 0.2, 0.1))
  par(omi=c(0.1, 0.3, 0.4, 0.1))
  
  layout(matrix(c(1:3, 4, 5:7, 4, 8:10, 4), nrow=4), heights=c(rep(1, 3), 0.3, 
                                                            rep(1, 3), 0.3, 
                                                            rep(1, 3), 0.3))
  #par(mfcol=c(3, 3))
  par(bty="n")
  
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    data_files_hist <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                         scale, "/historical/", vars[v]),
                                  pattern=paste0("Mean_", metrics[m]),
                                  full.names=TRUE, recursive=TRUE)
    #Read datasets
    data_files_rcp45 <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                          scale, "/rcp45/", vars[v]),
                                   pattern=paste0("Mean_", metrics[m]),
                                   full.names=TRUE, recursive=TRUE)
    
    data_files_rcp85 <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                          scale, "/rcp85/", vars[v]),
                                   pattern=paste0("Mean_", metrics[m]),
                                   full.names=TRUE, recursive=TRUE)
    
    
    
    #Should find 16 files in each case, check
    if (any(c(length(data_files_hist), length(data_files_rcp45),
              length(data_files_rcp85)) != 16)) {
      stop("Incorrect number of files found")
    }
    
    
    #Sources of uncertainty:
    #GCM, BC method and scenario
    
    
    #For each source, calculate absolute range for all possible combinations
    #Then take the average of these for plotting
    
    
    ### Full range ###
    
    #Calculate the full absolute range of future change
    
    #Load data
    data_hist  <- brick(lapply(data_files_hist, raster))
    data_rcp45 <- brick(lapply(data_files_rcp45, raster))
    data_rcp85 <- brick(lapply(data_files_rcp85, raster))
    
    #Calculate future change
    future_diff_rcp45 <- data_rcp45 - data_hist
    future_diff_rcp85 <- data_rcp85 - data_hist
    
    #Collate
    future_diff <- list(rcp45=future_diff_rcp45, rcp85=future_diff_rcp85)
    
    #Calculate full range (separately for each RCP, and then both RCPs)
    full_range_rcp45 <- max(future_diff_rcp45, na.rm=TRUE) - min(future_diff_rcp45, na.rm=TRUE)
    full_range_rcp85 <- max(future_diff_rcp85, na.rm=TRUE) - min(future_diff_rcp85, na.rm=TRUE)
    
    full_range_both_rcps <- max(brick(future_diff), na.rm=TRUE) - 
                            min(brick(future_diff), na.rm=TRUE)
    
    
    #Collate
    full_range <- list(rcp45=full_range_rcp45, rcp85=full_range_rcp85)
    
    rcps <- names(full_range)
    
    #Calculate areal mean of full range (for each RCP, and then both together)
    
    areal_mean_full_range <- list()
    for (r in rcps) areal_mean_full_range[[r]] <- areal_mean(full_range[[r]])

    areal_mean_full_range_both_rcps <- areal_mean(full_range_both_rcps)
    
    
    #Initialise
    plot_data <- list()
    
    plot_data_all <-list()
    
    aus_average <- list() #save continental average
    
    
    #Collate files by scenario
    gcm_files <- list(rcp45=data_files_rcp45, rcp85=data_files_rcp85)
    
    
    
    ###########
    ### GCM ###
    ###########
    
    plot_data[["GCM"]]     <- list()
    plot_data_all[["GCM"]] <- list()
    aus_average[["GCM"]]   <- list()
          
    #Four combinations for each RCP:
    # CCAM/ISIMIP/MRNBC/QME bc x 4 GCMS 
    
    bc_methods <- c("CCAM", "ISIMIP2b", "MRNBC", "QME")
    
    
    #Loop through scenarios
    for (g in rcps) {
      
      bc_data <- brick()
      
      #Loop through bc methods
      for (b in 1:length(bc_methods)) {
        
        #Find layers corresponding to BC method
        bc_ind <- which(grepl(bc_methods[b], gcm_files[[g]]))
        
        #Grab those
        all_data <- future_diff[[g]][[bc_ind]]
        
        #Calculate range for each bc layer
        bc_data <- addLayer(bc_data, max(all_data, na.rm=TRUE) - min(all_data, na.rm=TRUE))
        
        
      }

      #Save all the data
      plot_data_all$GCM[[g]] <- bc_data / full_range[[g]]
      
      
      #Save the mean 
      plot_data$GCM[[g]] <- mean(bc_data, na.rm=TRUE) / full_range[[g]]
      
      
      #Calculate continental averages (and save mean, min and max)
      area_mean <- areal_mean(bc_data) / areal_mean_full_range[[g]]
      
      aus_average$GCM[[g]] <- c(mean(area_mean), range(area_mean))
      
      
    } #RCPs
   
    
   
    
    ##################
    ### BC methods ###
    ##################
    
    plot_data[["BC"]]     <- list()
    plot_data_all[["BC"]] <- list()
    aus_average[["BC"]]   <- list()
    
    #Four combinations for each BC method
    #  CNRM/ACCESS/MIROC/GFLD gcm x 4 BC methods 
    
    gcms <- c("CNRM-CERFACS-CNRM-CM5", "CSIRO-BOM-ACCESS1-0",
              "MIROC-MIROC5", "NOAA-GFDL-GFDL-ESM2M")
    
    
    #Loop through scenarios
    for (g in rcps) {
      
      gcm_data <- brick()
      
      #Loop through GCMs
      for (b in 1:length(gcms)) {
        
        #Find layers corresponding to BC method
        gcm_ind <- which(grepl(gcms[b], gcm_files[[g]]))
        
        #Grab those
        all_data <- future_diff[[g]][[gcm_ind]]
        
        #Calculate range for each bc layer
        gcm_data <- addLayer(gcm_data, max(all_data, na.rm=TRUE) - min(all_data, na.rm=TRUE))
        
        
      }
      
      #Save all the data
      plot_data_all$BC[[g]] <- gcm_data / full_range[[g]]
      
      
      #Save the mean 
      plot_data$BC[[g]] <- mean(gcm_data) / full_range[[g]]
      
      
      #Calculate continental averages (and save mean, min and max)
      area_mean <- areal_mean(gcm_data) / areal_mean_full_range[[g]]
      
      aus_average$BC[[g]] <- c(mean(area_mean), range(area_mean))
      
      
    } #RCPs
    
    
    
    
    
    #####################
    ### RCP scenarios ###
    #####################
    
    
    plot_data[["RCP"]]     <- list()
    plot_data_all[["RCP"]] <- list()
    aus_average[["RCP"]]   <- list()
    
    #16 combinations
    #(each bc method + each GCM x 2 RCP scenarios)
    
    
    rcp_data <- brick()
    
    #Loop through scenarios
    for (b in bc_methods) {
      
      #Loop through GCMs
      for (g in gcms) {
   
        #Find layers corresponding to BC method
        rcp_ind <- lapply(gcm_files, function(f) which(grepl(paste0(b, "_", g), f)))
        
        #Grab those
        all_data <- brick(lapply(1:length(rcp_ind), function(x) future_diff[[x]][[rcp_ind[[x]]]]))
        
        #Calculate range for each bc layer
        rcp_data <- addLayer(rcp_data, max(all_data, na.rm=TRUE) - min(all_data, na.rm=TRUE))
        
    
      }
    }
    
    
    
    #Save all the data
    plot_data_all$RCP <- rcp_data / full_range_both_rcps
    
    
    #Save the mean 
    plot_data$RCP <- mean(rcp_data) / full_range_both_rcps
    
    
    #Calculate continental averages (and save mean, min and max)
    area_mean <- areal_mean(rcp_data) / areal_mean_full_range_both_rcps 

    aus_average$RCP <- c(mean(area_mean), range(area_mean))
    
    
    
    
    
    ################
    ### Plotting ###
    ################
    
    
    #Only plot RCP4.5 here for GCM and BC. Making a separate SI figure for RCP8.5
    
    plot_col <- cols(length(lims)-1)
    
    
    ### GCM ###
    
    #Plot (multiply by 100 to convert to percentage)
    image(plot_data$GCM$rcp45 * 100, breaks=lims, 
          col=plot_col, axes=FALSE, ann=FALSE, asp=1)
    
    
    #Australia outline
    map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
    
    
    #Variable label
    mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
    
    
    #Source of uncertainty label   
    if (v==1) mtext(side=2, line=1, text="GCM")
    
  
    #Aus average
    text <- paste0(round(aus_average$GCM$rcp45[1]*100), "%")
    
    ci_text <- paste0("[", round(aus_average$GCM$rcp45[2] * 100), "-",
                      round(aus_average$GCM$rcp45[3] * 100), "%]")
    
    
    text(118, -12, text, cex=0.8)
    text(118, -15, ci_text, cex=0.8, adj=0.5)

    # mtext(side=3, line=-1, text=text, adj=0.0, cex=0.45)
    # mtext(side=3, line=-3, text=text1, adj=-0.05, cex=0.45)
    # 
    
    
    
    ### BC methods ###
    
    #Plot (multiply by 100 to convert to percentage)
    image(plot_data$BC$rcp45 * 100, breaks=lims, 
          col=plot_col, axes=FALSE, ann=FALSE, asp=1)
    
    
    #Australia outline
    map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
    
    
    #Source of uncertainty label   
    if (v==1) mtext(side=2, line=1, text="BC method")
    
    #Aus average
    text <- paste0(round(aus_average$BC$rcp45[1]*100), "%")
    
    ci_text <- paste0("[", round(aus_average$BC$rcp45[2] * 100), "-",
                      round(aus_average$BC$rcp45[3] * 100), "%]")
    
    
    text(118, -12, text, cex=0.8)
    text(118, -15, ci_text, cex=0.8, adj=0.5)
    
    
    
    
    ### RCP scenario ###
    
    #Plot (multiply by 100 to convert to percentage)
    image(plot_data$RCP * 100, breaks=lims, 
          col=plot_col, axes=FALSE, ann=FALSE, asp=1)
    
    
    #Australia outline
    map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
    
    
    #Source of uncertainty label   
    if (v==1) mtext(side=2, line=1, text="RCP scenario")
    
    
    #Aus average
    text <- paste0(round(aus_average$RCP[1]*100), "%")
    
    ci_text <- paste0("[", round(aus_average$RCP[2] * 100), "-",
                      round(aus_average$RCP[3] * 100), "%]")
    
    
    text(118, -12, text, cex=0.8)
    text(118, -15, ci_text, cex=0.8, adj=0.5)
    
    
    
    
    #Legend 
    if (v == 1) {
      
      #Empty plot
      plot(1, type="n", yaxt="n", xaxt="n")
      
      #Legend
      len <- length(lims)-1
      add_raster_legend2(cols=plot_col, limits=lims[2:len],
                         main_title="%", plot_loc=c(0.35,0.65,0.63, 0.77), 
                         title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
      
      
    }
      
    
  } #variables
  
  dev.off()
} #metrics







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
vars <- c("pr", "qtot", "sm_root")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("frequency", "duration", "rel_intensity") #, "frequency")


#Experiments
exp <- c("rcp45", "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")


outdir <- paste0(path, "/Figures")
dir.create(outdir)




#Loop through metrics
for (m in 1:length(metrics)) {
  

  
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

    #Calculate future change
    future_diff_rcp45 <- data_rcp45 - data_hist

    #Collate
    future_diff <- list(rcp45=future_diff_rcp45)
    
    #Calculate full range (separately for each RCP, and then both RCPs)
    full_range_rcp45 <- max(future_diff_rcp45, na.rm=TRUE) - min(future_diff_rcp45, na.rm=TRUE)

    full_range_both_rcps <- max(brick(future_diff), na.rm=TRUE) - 
      min(brick(future_diff), na.rm=TRUE)
    
    
    #Collate
    full_range <- list(rcp45=full_range_rcp45)
    
    rcps <- names(full_range)
    
    #Calculate areal mean of full range (for each RCP, and then both together)
    
    areal_mean_full_range <- list()
    for (r in rcps) areal_mean_full_range[[r]] <- areal_mean(full_range[[r]])
    
    areal_mean_full_range_both_rcps <- areal_mean(full_range_both_rcps)
    
    
    #Initialise
    plot_data <- list()
    
    plot_data_all <-list()
    
    aus_average <- list() #save continental average
    aus_average_no_ccam <- list() #save continental average
    
    
    #Collate files by scenario
    gcm_files <- list(rcp45=data_files_rcp45)
    
    
    
    ###########
    ### GCM ###
    ###########
    
    plot_data[["GCM"]]     <- list()
    plot_data_all[["GCM"]] <- list()
    aus_average[["GCM"]]   <- list()
    
    #Four combinations for each RCP:
    # CCAM/ISIMIP/MRNBC/QME bc x 4 GCMS 
    
    bc_methods <- c("CCAM", "ISIMIP2b", "MRNBC", "QME")
    
    
 
    
    ##################
    ### BC methods ###
    ##################
    
    plot_data[["BC"]]     <- list()
    plot_data_all[["BC"]] <- list()
    aus_average[["BC"]]   <- list()
    aus_average_no_ccam[["BC"]]   <- list()
    
    #Four combinations for each BC method
    #  CNRM/ACCESS/MIROC/GFLD gcm x 4 BC methods 
    
    gcms <- c("CNRM-CERFACS-CNRM-CM5", "CSIRO-BOM-ACCESS1-0",
              "MIROC-MIROC5", "NOAA-GFDL-GFDL-ESM2M")
    
    
    #Loop through scenarios
    for (g in rcps) {
      
      gcm_data         <- brick()
      gcm_data_no_ccam <- brick()
      
      #Loop through GCMs
      for (b in 1:length(gcms)) {
        
        ### ALL ###
        
        #Find layers corresponding to BC method
        gcm_ind <- which(grepl(gcms[b], gcm_files[[g]]))
        
        #Grab those
        all_data <- future_diff[[g]][[gcm_ind]]
        
        #Calculate range for each bc layer
        gcm_data <- addLayer(gcm_data, max(all_data, na.rm=TRUE) - min(all_data, na.rm=TRUE))
        
        ### NO CCAM ###
        
        ind_no_ccam <- gcm_ind[!grepl("/CCAM/", gcm_files$rcp45[gcm_ind])]
        
        #Grab those
        all_data_no_ccam <- future_diff[[g]][[ind_no_ccam]]
        
        #Calculate range for each bc layer
        gcm_data_no_ccam <- addLayer(gcm_data_no_ccam, max(all_data_no_ccam, na.rm=TRUE) - 
                                       min(all_data_no_ccam, na.rm=TRUE))
        
        
        
        
        
      }
      
      #Save all the data
      #plot_data_all$BC[[g]] <- gcm_data / full_range[[g]]
      
      
      #Save the mean 
      #plot_data$BC[[g]] <- mean(gcm_data) / full_range[[g]]
      
      
      #Calculate continental averages (and save mean, min and max)
      area_mean <- areal_mean(gcm_data) / areal_mean_full_range[[g]]
      
      aus_average$BC[[g]] <- c(mean(area_mean), range(area_mean))
      
      area_mean_no_ccam <- areal_mean(gcm_data_no_ccam) / areal_mean_full_range[[g]]
      aus_average_no_ccam$BC_no_ccam[[g]] <- c(mean(area_mean_no_ccam), range(area_mean_no_ccam))
      
      
      
    } #RCPs
    
    
    print(paste0("Metric: ", metrics[m], ", Variable: ", vars[v]))
    print(paste0("Mean all: ", round(mean(area_mean), 2)))
    print(paste0("Mean no CCAM: ", round(mean(area_mean_no_ccam), 2)))
    print(paste0("difference: ", round(mean(area_mean - area_mean_no_ccam), 2)))
    
    
    
   
 
    #Aus average
    text <- paste0(round(aus_average$BC$rcp45[1]*100), "%")
    
    ci_text <- paste0("ALL methods: [", round(aus_average$BC$rcp45[2] * 100), "-",
                      round(aus_average$BC$rcp45[3] * 100), "%]")
    
    #print(ci_text)
  
     } #variables
  
} #metrics







library(raster)
library(RColorBrewer)
library(maptools)
library(maps)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/"


#Source functions
source(paste0(path,"/scripts/R/functions/wtd_stdev.R"))
source(paste0(path,"/scripts/R/functions/add_raster_legend.R"))



#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3

#baseline   <- "1970-2005_vs_2064-2099"


#Variables
vars <- c("qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Runoff", "Soil moisture") #labels for plotting

gcm_labels <- list("CNRM-CERFACS-CNRM-CM5" = "CNRM-CM5",
                   "CSIRO-BOM-ACCESS1-0"   = "ACCESS1-0",
                   "MIROC-MIROC5"          = "MIROC5",
                   "NOAA-GFDL-GFDL-ESM2M"  = "GFDL-ESM2M")



#List metrics
metrics <- c("duration", "rel_intensity", "timing")#, "frequency")


#Raster 1x1 deg for resampling to common resolution
#extent_raster  <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)

outdir <- paste0(path, "/Figures_CABLE_AWRA/")
dir.create(outdir)


#####################
### Plot settings ###
#####################


#Set plot colours

#Historical mean
cols_hist <- colorRampPalette(c("#ffffe5", "#fee391",
                                         "#fe9929", "#cc4c02"))
                                         

#Future difference
col_opts <- rev(brewer.pal(11, 'RdYlBu'))
#col_opts[6] <- "grey80"
cols_diff <- colorRampPalette(col_opts) 

#Limits
lims_hist <- list(pr =   list(duration  = c(1, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 1000),
                              rel_intensity = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                              frequency = c(0, 8.5, 8.8, 9.1, 9.4, 9.7, 10, 10.3, 1000)),
                  qtot = list(duration  = c(2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 1000),
                              rel_intensity = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                              frequency = c(0, 8.5, 8.8, 9.1, 9.4, 9.7, 10, 10.3, 1000)),
                  sm =   list(duration  = c(2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 1000),
                              rel_intensity = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                              frequency = c(0, 8.5, 8.8, 9.1, 9.4, 9.7, 10, 10.3, 1000)))

lims_diff <- list(duration      = c(-1000, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 1000),
                  intensity     = c(-10000, -30, -25, -20, -16, -12, -8, -4, -2, 0, 2, 4, 8, 12, 16, 20, 25, 30, 10000),
                  rel_intensity = c(-10000, -8, -6, -4, -2, 0,  2, 4, 6, 8, 10000), #c(-100, -40, -30, -20, -10, 10, 20, 30, 40, 100),
                  frequency     = c(-100, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5,0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 100),
                  timing        = c(-10000, -20, -15, -10, -5, 0,  5, 10, 15, 20, 10000))


unit <- c("months", "% points", "% points") # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_CABLE_AWRA_changes_by_GCM_in_", metrics[m], "_",
               percentile, "_", scale, "_", vars[v], ".png"),
        height=10.5, width=8.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.3, 0.4, 0.1))
    
    layout(matrix(c(1:12, 13, 13, 13), nrow=5, byrow=TRUE), heights=c(1, 1, 1, 1, 0.3,
                                                                  1 ,1, 1, 1, 0.3,
                                                                  1 ,1, 1, 1, 0.3))
    
    
    #par(mfcol=c(3, 3))
    par(bty="n")
    
    
    
    ### Read datasets ###
    
    
    #List GCMs
    gcms <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                              scale, "/", vars[v], "/CO2/"))
    
    
    for (g in 1:length(gcms)) {
      
      
        
      
      #############
      ### CABLE ###
      #############
      
      
      if (metrics[m] == "timing") {
        
        ############
        ### AWRA ###
        ############
        
        
        data_files_awra <- list.files(paste0(path, "/drought_metrics/", 
                                           scale, "-month/MRNBC/", gcms[g]),
                                           pattern=glob2rx(paste0("drought_metrics_*", vars[v], "*rcp45*.nc")),
                                           full.names=TRUE, recursive=TRUE)
    
        #CO2
        data_files_CO2  <- list.files(paste0(path, "/drought_metrics_CABLE/", 
                                                  scale, "-month/CO2/", gcms[g]),
                                           pattern=glob2rx(paste0("drought_metrics_CABLE_CO2_*", vars[v], 
                                                                  ".nc")),
                                           full.names=TRUE)
        
        
        #noCO2
        data_files_noCO2  <- list.files(paste0(path, "/drought_metrics_CABLE/", 
                                             scale, "-month/noCO2/", gcms[g]),
                                      pattern=glob2rx(paste0("drought_metrics_CABLE_noCO2_*", vars[v], 
                                                             ".nc")),
                                      full.names=TRUE)
        
        
        #Check number of files
        if (any (c(length(data_files_noCO2), length(data_files_CO2), length(data_files_awra)) != 1)){
          stop("Wrong number of files")
        }
        
        
        #Load data
        data_awra   <- brick(data_files_awra, varname=metrics[m])
        data_CO2    <- brick(data_files_CO2, varname=metrics[m])
        data_noCO2  <- brick(data_files_noCO2, varname=metrics[m])
        
        #Calculate historical and future mean
        
        #First get date indices so that use 1970-2005 for historical
        #and 2064-2099 for future
        
        dates <- getZ(data_awra)
        
        start_hist <- which(dates  == "1970-01-16")
        end_hist <- which(dates  == "2005-12-16")
        
        len_time <- length(start_hist:end_hist)
        
        start_rcp45 <- which(dates  == "2064-01-16")
        end_rcp45 <- which(dates  == "2099-12-16")
        
        #Historical mean (percentage of time in drought)
        data_hist_awra  <- sum(data_awra[[start_hist:end_hist]], na.rm=TRUE) / len_time * 100
        data_hist_CO2   <- sum(data_CO2[[start_hist:end_hist]], na.rm=TRUE) / len_time *100
        data_hist_noCO2 <- sum(data_noCO2[[start_hist:end_hist]], na.rm=TRUE) / len_time *100
        
        #Future mean
        data_rcp45_awra  <- sum(data_awra[[start_rcp45:end_rcp45]], na.rm=TRUE) / len_time * 100
        data_rcp45_CO2   <- sum(data_CO2[[start_rcp45:end_rcp45]], na.rm=TRUE) / len_time *100
        data_rcp45_noCO2 <- sum(data_noCO2[[start_rcp45:end_rcp45]], na.rm=TRUE) / len_time *100
      
        
      } else {
        
        
        ############
        ### AWRA ###
        ############
        
        
        data_files_hist_awra <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                                  scale, "/historical/", vars[v], "/MRNBC/", gcms[g]),
                                           pattern=paste0("Mean_", metrics[m]),
                                           full.names=TRUE, recursive=TRUE)
        
        #Read datasets
        data_files_rcp45_awra <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                                   scale, "/rcp45/", vars[v],"/MRNBC/", gcms[g]),
                                            pattern=paste0("Mean_", metrics[m]),
                                            full.names=TRUE, recursive=TRUE)
        
        #Check number of files
        if (any (c(length(data_files_hist_awra), length(data_files_rcp45_awra)) != 1)){
          stop("Wrong number of CO2 files")
        }
        
        #CO2
        data_files_hist_CO2  <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                  scale, "/", vars[v], "/CO2/", gcms[g]),
                                           pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                  "*_CO2_*_historical_*.nc")),
                                           full.names=TRUE, recursive=TRUE)
        
        data_files_rcp45_CO2 <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                  scale, "/", vars[v], "/CO2/", gcms[g]),
                                           pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                  "*_CO2_*_rcp45_*.nc")),
                                           full.names=TRUE, recursive=TRUE)
        
        #Check number of files
        if (any (c(length(data_files_hist_CO2), length(data_files_rcp45_CO2)) != 1)){
          stop("Wrong number of CO2 files")
        }
        
        
        #noCO2
        data_files_hist_noCO2  <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                    scale, "/", vars[v], "/noCO2/", gcms[g]),
                                             pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                    "*_noCO2_*_historical_*.nc")),
                                             full.names=TRUE, recursive=TRUE)
        
        data_files_rcp45_noCO2 <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                    scale, "/", vars[v], "/noCO2/", gcms[g]),
                                             pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                    "*_noCO2_*_rcp45_*.nc")),
                                             full.names=TRUE, recursive=TRUE)
        
        #Check number of files
        if (any (c(length(data_files_hist_noCO2), length(data_files_rcp45_noCO2)) != 1)){
          stop("Wrong number of CO2 files")
        }
        
        
        #Load data
        data_hist_awra    <- raster(data_files_hist_awra)
        data_hist_CO2     <- raster(data_files_hist_CO2)
        data_hist_noCO2   <- raster(data_files_hist_noCO2)
        
        data_rcp45_awra   <- raster(data_files_rcp45_awra)
        data_rcp45_CO2    <- raster(data_files_rcp45_CO2)
        data_rcp45_noCO2  <- raster(data_files_rcp45_noCO2)
        
      }

      
      
      #Calculate future change
      future_diff_rcp45_awra   <- data_rcp45_awra - data_hist_awra
      future_diff_rcp45_CO2    <- data_rcp45_CO2 - data_hist_CO2
      future_diff_rcp45_noCO2  <- data_rcp45_noCO2 - data_hist_noCO2
      
      
      
      ################
      ### Plotting ###
      ################
      
    
      ### Future change ###
      
      plot_data <- list(awra  = future_diff_rcp45_awra,
                        CO2   = future_diff_rcp45_CO2,
                        noCO2 = future_diff_rcp45_noCO2)
      
   
      labs <- c("AWRA", "CABLE CO2", "CABLE noCO2")   
      
      
      #Loop through experiments
      for (p in 1:length(plot_data)) {
        
        
        #Plot
        len <- length(lims_diff[[metrics[m]]])
        image(plot_data[[p]], breaks=lims_diff[[metrics[m]]], 
              col=c(cols_diff(len-1)),
              axes=FALSE, ann=FALSE, asp=1)
        
        
        #Australia outline
        map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
        
        
        #Model label
        if(g==1) mtext(side=3, line=1, font=2, text=labs[p], xpd=NA)
        
        
        if (p==1) mtext(side=2, line=1, text=gcm_labels[[gcms[g]]])
        
        
      }
      
      
   
      #Empty plot
      if (g == length(gcms)) {
        
        plot(1, type="n", bty="n", yaxt="n", xaxt="n")
      
      #Legend
      len1 <- length(lims_diff[[metrics[m]]])-1
      add_raster_legend2(cols=cols_diff(len1), limits=lims_diff[[metrics[m]]][2:len1],
                         main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                         title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
      }
      
      
    } #gcms
    
    dev.off ()
    
  } #metrics
  
} #variables







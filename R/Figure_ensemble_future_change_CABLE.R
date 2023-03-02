library(raster)
library(RColorBrewer)
library(maptools)
library(maps)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"


#Source functions
source(paste0(path,"/scripts/R/functions/perc_agree_on_sign_weighted_AR4_method.R"))
source(paste0(path,"/scripts/R/functions/perc_agree_on_mean_weighted.R"))
source(paste0(path,"/scripts/R/functions/perc_agree_with_obs_mean_weighted.R"))
source(paste0(path,"/scripts/R/functions/wtd_stdev.R"))
source(paste0(path,"/scripts/R/functions/add_raster_legend.R"))



#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3

#baseline   <- "1970-2005_vs_2064-2099"


#Variables
vars <- c("qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")


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
                  frequency     = c(-100, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5,0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 100))


unit <- c("months", "% points") # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  ### Set up figure ###
  png(paste0(outdir, "/FigureX", "_CABLE_Mean_changes_in_", metrics[m], "_",
             percentile, "_", scale, ".png"),
      height=10.5, width=6.3, units="in", res=400)
  
  
  par(mai=c(0, 0.1, 0.2, 0.1))
  par(omi=c(0.1, 0.3, 0.4, 0.1))
  
  layout(matrix(c(1:8, 3, 9:10, 6), nrow=6), heights=c(1 ,1, 0.3, 1, 1, 0.3,
                                                       1 ,1, 0.3, 1, 1, 0.3))
                                                       

  #par(mfcol=c(3, 3))
  par(bty="n")
  
  

  #Loop through variables
  for (v in 1:length(vars)) {
  
    
    ### Read datasets ###
    
    #CO2
    data_files_hist_CO2  <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                          scale, "/", vars[v], "/CO2/"),
                                       pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                              "*_CO2_*_historical_*.nc")),
                                  full.names=TRUE, recursive=TRUE)
    
    
    data_files_rcp45_CO2 <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                          scale, "/", vars[v], "/CO2/"),
                                       pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                              "*_CO2_*_rcp45_*.nc")),
                                   full.names=TRUE, recursive=TRUE)

    #Check number of files
    if (any (c(length(data_files_hist_CO2), length(data_files_rcp45_CO2)) != 4)){
      stop("Wrong number of CO2 files")
    }
    
    
    #noCO2
    data_files_hist_noCO2  <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                scale, "/", vars[v], "/noCO2/"),
                                         pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                "*_noCO2_*_historical_*.nc")),
                                         full.names=TRUE, recursive=TRUE)
    
    data_files_rcp45_noCO2 <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                              scale, "/", vars[v], "/noCO2/"),
                                       pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                              "*_noCO2_*_rcp45_*.nc")),
                                       full.names=TRUE, recursive=TRUE)
    
    #Check number of files
    if (any (c(length(data_files_hist_noCO2), length(data_files_rcp45_noCO2)) != 4)){
      stop("Wrong number of CO2 files")
    }
    
    
    #Load data
    data_hist_CO2    <- brick(lapply(data_files_hist_CO2, raster))
    data_hist_noCO2  <- brick(lapply(data_files_hist_noCO2, raster))
    
    data_rcp45_CO2   <- brick(lapply(data_files_rcp45_CO2, raster))
    data_rcp45_noCO2 <- brick(lapply(data_files_rcp45_noCO2, raster))
    
    
    
    #Calculate historical model mean
    ens_median_hist_CO2   <- calc(data_hist_CO2, median)
    ens_median_hist_noCO2 <- calc(data_hist_noCO2, median)
    
    
    #Calculate future change
    future_diff_rcp45_CO2   <- data_rcp45_CO2 - data_hist_CO2
    future_diff_rcp45_noCO2 <- data_rcp45_noCO2 - data_hist_noCO2
    
    
    #Ensemble mean of future change
    ens_median_rcp45_CO2   <- calc(future_diff_rcp45_CO2, median)
    ens_median_rcp45_noCO2 <- calc(future_diff_rcp45_noCO2, median)
    

    #Calculate model agreement
    mod_agr_rcp45_CO2 <- calc(future_diff_rcp45_CO2, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                      weights=rep(1, nlayers(data_rcp45_CO2))))
    
    mod_agr_rcp45_noCO2 <- calc(future_diff_rcp45_noCO2, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                                     weights=rep(1, nlayers(data_rcp45_noCO2))))
    
    
      
  
    ################
    ### Plotting ###
    ################
    
    
    ### Historical median ###
    
    plot_col <- cols_hist(length(lims_hist[[vars[v]]][[metrics[m]]])-1)
    
    
    ### Rising CO2 ###
    
    #Plot
    image(ens_median_hist_CO2, breaks=lims_hist[[vars[v]]][[metrics[m]]], 
          col=plot_col, axes=FALSE, ann=FALSE, asp=1)
      
    #Australia outline
    map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
    
    
    #Variable label
    mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
    
    
    #Experiment label   
    if (v==1) mtext(side=2, line=1, text="CO2 historical")
    
    
    ### constant CO2 ###
    
    #Plot
    image(ens_median_hist_noCO2, breaks=lims_hist[[vars[v]]][[metrics[m]]], 
          col=plot_col, axes=FALSE, ann=FALSE, asp=1)
    
    #Australia outline
    map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
    
    #Experiment label   
    if (v==1) mtext(side=2, line=1, text="noCO2 historical")
    
    
    ### Legend ###
    
    
    #Legend
    if (v == 1) {
      
      #Empty plot
      plot(1, type="n", yaxt="n", xaxt="n")
      
      len <- length(lims_hist[[vars[v]]][[metrics[m]]])-1
      add_raster_legend2(cols=plot_col, limits=lims_hist[[vars[v]]][[metrics[m]]][2:len],
                         main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                         title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
    }
     
    
    
    ### Future change ###
    
    plot_data <- list(CO2 = ens_median_rcp45_CO2,
                      noCO2 = ens_median_rcp45_noCO2)
    
    agr_data <- list(CO2 = mod_agr_rcp45_CO2,
                     noCO2 = mod_agr_rcp45_noCO2)
    
    co2_labs <- c("CO2 RCP-4.5", "noCO2 RCP-4.5")   
    
    
    #Loop through experiments
    for (p in 1:length(plot_data)) {
      
      
      #Set pixels with no model agreement to some crazy values
      #(will plot these in grey)
      mask_value <- -1000000000

      #plot_data[[p]][agr_data[[p]] == 0] <- mask_value

      
      #Plot
      len <- length(lims_diff[[metrics[m]]])
      image(plot_data[[p]], breaks=c(mask_value, mask_value+1, lims_diff[[metrics[m]]][2:len]), 
            col=c("grey80", cols_diff(len-1)),
            axes=FALSE, ann=FALSE, asp=1)
      
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"

      
      #Variable label
      #if(p==1) mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
    
     
      #Calculate and print land area where model changes are robust
      area_agree <- area(agr_data[[p]])
      area_agree[agr_data[[p]] != 1 | is.na(agr_data[[p]])] <- NA
      
      area_total <- mask(area(agr_data[[p]]), agr_data[[p]])
      
      land_area <- (sum(values(area_agree), na.rm=TRUE)  / 
                   sum(values(area_total), na.rm=TRUE)) *100
      
      
      mtext(side=3, line=-3, text=paste0(sprintf("%.1f", land_area), "%"), adj=0.1, cex=0.65)
      
      
      if (v==1) mtext(side=2, line=1, text=co2_labs[p])
      
      
    }
    
    
    if (v==1) {
      
      #Empty plot
      plot(1, type="n", bty="n", yaxt="n", xaxt="n")
      
      #Legend
      len <- length(lims_diff[[metrics[m]]])-1
      add_raster_legend2(cols=cols_diff(len), limits=lims_diff[[metrics[m]]][2:len],
                                         main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                                         title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
    }
     
    
  } #metrics
  
  dev.off ()
} #variables







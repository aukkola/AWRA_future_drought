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


#Load world shapefile
world <- readShapePoly(paste0(path, "../World_shapefile/World"))


#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3

#baseline   <- "1970-2005_vs_2064-2099"


#Variables
vars <- c("pr", "qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity", "frequency")


#Experiments
exp <- c("rcp45", "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Raster 1x1 deg for resampling to common resolution
#extent_raster  <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)

outdir <- paste0(path, "/Figures")
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
                  frequency     = c(-100, -20, -15, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 15, 20, 100))


unit <- c("months", "% points", "% points") # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  if (metrics[m] == "frequency") {
    height= 4
  } else {
    height=6.5
  }
  
  ### Set up figure ###
  png(paste0(outdir, "/FigureX", "_Mean_changes_in_", metrics[m], "_",
             percentile, "_", scale, ".png"),
      height=height, width=7.3, units="in", res=400)
  
  
  par(mai=c(0, 0.05, 0.2, 0.0))
  par(omi=c(0.1, 0.3, 0.4, 0.1))
  
   #par(mfcol=c(3, 3))
  par(bty="n")
  
  
  if (metrics[m] == "frequency") {
    layout(matrix(c(1:3, 4:5, 3, 6:7, 3), nrow=3), heights=c(1 , 1, 0.3, 
                                                             1, 1, 0.3,
                                                             1, 1, 0.3))
  } else {
    layout(matrix(c(1:4, 5, 6:9, 5, 10:13, 5), nrow=5), heights=c(1 ,0.3, 1, 1, 0.3,
                                                                  1 ,0.3, 1, 1, 0.3,
                                                                  1 ,0.3, 1, 1, 0.3))
  }
  

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
    
    
    data_hist  <- brick(lapply(data_files_hist, raster))
    data_rcp45 <- brick(lapply(data_files_rcp45, raster))
    data_rcp85 <- brick(lapply(data_files_rcp85, raster))
    
    
    #Calculate historical model mean
    ens_median_hist <- calc(data_hist, median)
    
    
    #Calculate future change
    
    future_diff_rcp45 <- data_rcp45 - data_hist
    future_diff_rcp85 <- data_rcp85 - data_hist
    
    #Ensemble mean of future change
    ens_median_rcp45 <- calc(future_diff_rcp45, median)
    ens_median_rcp85 <- calc(future_diff_rcp85, median)
    

    #Calculate model agreement
    mod_agr_rcp45 <- calc(future_diff_rcp45, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                      weights=rep(1, nlayers(data_rcp45))))
    mod_agr_rcp85 <- calc(future_diff_rcp85, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                      weights=rep(1, nlayers(data_rcp85))))
    
      
  
    ################
    ### Plotting ###
    ################
    
    if (metrics[m] != "frequency") {
      
      ### Historical median ###
      
      plot_col <- cols_hist(length(lims_hist[[vars[v]]][[metrics[m]]])-1)
      
      #Plot
      image(ens_median_hist, breaks=lims_hist[[vars[v]]][[metrics[m]]], 
            col=plot_col, axes=FALSE, ann=FALSE, asp=1)
      
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
      
      
      #Variable label
      mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
      
      
      #Experiment label   
      if (v==1) mtext(side=2, line=1, text="Historical")
      
      
      #Empty plot
      plot(1, type="n", yaxt="n", xaxt="n")
      
      #Legend
      len <- length(lims_hist[[vars[v]]][[metrics[m]]])-1
      add_raster_legend2(cols=plot_col, limits=lims_hist[[vars[v]]][[metrics[m]]][2:len],
                         main_title=unit[m], plot_loc=c(0.1,0.9,0.63, 0.77), 
                         title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
      
      
      
    }
    
    ### Future change ###
    
    plot_data <- list(rcp45 = ens_median_rcp45,
                      rcp85 = ens_median_rcp85)
    
    agr_data <- list(rcp45 = mod_agr_rcp45,
                     rcp85 = mod_agr_rcp85)
    
    
    #Loop through experiments
    for (p in 1:length(plot_data)) {
      
      
      
      
      #Set pixels with no model agreement to some crazy values
      #(will plot these in grey)
      mask_value <- -1000000000
        
      plot_data[[p]][agr_data[[p]] == 0] <- mask_value
      
      
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
      
      
      # #Stippling (where models don't agree)
      # 
      # stipple_raster <- raster(resolution=c(0.25, 0.25), ext=extent(agr_data[[p]]))
      # 
      # #Resample so points larger
      # fut_mod_agr <- resample(agr_data[[p]], stipple_raster)
      # 
      # 
      # #Find pixels where model's don't agree
      # ind    <- which(values(fut_mod_agr) == 0)
      # coords <- coordinates(fut_mod_agr)
      # 
      # #Add stippling
      # points(coords[ind,1], coords[ind,2], pch=20, lwd=lwd, cex=cex)
      # 
      # #Panel number
      # if (p ==1 ){
      #   mtext(side=3, line=0, adj=0, text="c", font=2, xpd=NA, cex=panel_cex) #panel number
      #   
      # } else {
      #   mtext(side=3, line=0, adj=0, text="d", font=2, xpd=NA, cex=panel_cex) #panel number
      #   
      # }
      
      if (v==1) mtext(side=2, line=1, text=exp_labels[p])
      
      
      #Variable label
      if (metrics[m] == "frequency" & p==1) mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
      
        
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







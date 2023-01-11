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

scale      <- "scale_3"

baseline   <- "1970-2005_vs_2064-2099"


#Variables
vars <- c("pr", "qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")


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
col_opts[6] <- "grey80"
cols_diff <- colorRampPalette(col_opts) 

#Limits
lims_hist <- list(duration  = c(1, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 1000),
                  intensity = c(0, 10, 15, 20, 25, 30, 35, 40, 10000),
                  rel_intensity = c(0, 10, 15, 20, 25, 30, 35, 40, 100),
                  frequency = c(0, 8.5, 8.8, 9.1, 9.4, 9.7, 10, 10.3, 1000))

lims_diff <- list(duration      = c(-1000, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 1000),
                  intensity     = c(-10000, -30, -25, -20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20, 25, 30, 10000),
                  rel_intensity = c(-100, -40, -30, -20, -10, 10, 20, 30, 40, 100),
                  frequency     = c(-100, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 100))


unit <- c("months", expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  ### Set up figure ###
  png(paste0(outdir, "/FigureX", "_Mean_changes_in_", metrics[m], "_",
             percentile, "_", scale, ".png"),
      height=5.5, width=8.3, units="in", res=400)
  
  
  par(mai=c(0, 0.1, 0, 0.1))
  par(omi=c(0.6, 0.3, 0.5, 0.1))
  
  par(mfcol=c(2, 3))
  par(bty="n")
  
  

  #Loop through variables
  for (v in 1:length(vars)) {
  
     
    #Read datasets
    data_files_rcp45 <- list.files(paste0(path, "/Mean_drought_metrics/rcp45/",
                                          vars[v], "/", metrics[m]),
                                   pattern=paste0(scale, "_", baseline),
                                   full.names=TRUE)

    data_files_rcp85 <- list.files(paste0(path, "/Mean_drought_metrics/rcp85/",
                                          vars[v], "/", metrics[m]),
                                   pattern=paste0(scale, "_", baseline),
                                   full.names=TRUE)
    
    
    data_rcp45 <- brick(lapply(data_files_rcp45, raster))
    data_rcp85 <- brick(lapply(data_files_rcp85, raster))
    
    
    #Calculate ensemble median
    ens_median_rcp45 <- calc(data_rcp45, median)
    ens_median_rcp85 <- calc(data_rcp85, median)
    

    #Calculate model agreement
    mod_agr_rcp45 <- calc(data_rcp45, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                      weights=rep(1, nlayers(data_rcp45))))
    mod_agr_rcp85 <- calc(data_rcp85, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                      weights=rep(1, nlayers(data_rcp85))))
    
      
  
    ################
    ### Plotting ###
    ################
  
    plot_data <- list(rcp45 = ens_median_rcp45,
                      rcp85 = ens_median_rcp85)
    
    agr_data <- list(rcp45 = mod_agr_rcp45,
                     rcp85 = mod_agr_rcp85)
    
    
    #Loop through experiments
    for (p in 1:length(plot_data)) {
      
      
      
      #Plot
      image(plot_data[[p]], breaks=lims_diff[[metrics[m]]], col=cols_diff(length(lims_diff[[metrics[m]]])-1),
            axes=FALSE, ann=FALSE, asp=1)
      
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"

      
      #Variable label
      if(p==1) mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
      
      
      #Legend
      len <- length(lims_diff[[metrics[m]]])-1
      if(p==2 & v==2) add_raster_legend2(cols=cols_diff(len), limits=lims_diff[[metrics[m]]][2:len],
                                  main_title=unit[m], plot_loc=c(0.12,0.88,-0.17,-0.14), 
                                  title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
      
     
      #Calculate and print land area where model changes are robust
      area_agree <- area(agr_data[[p]])
      area_agree[agr_data[[p]] != 1 | is.na(agr_data[[p]])] <- NA
      
      area_total <- mask(area(agr_data[[p]]), agr_data[[p]])
      
      land_area <- (sum(values(area_agree), na.rm=TRUE)  / 
                   sum(values(area_total), na.rm=TRUE)) *100
      
      
      mtext(side=3, line=-5, text=paste0(sprintf("%.1f", land_area), "%"), adj=0.1, cex=0.65)
      
      
      #Stippling (where models don't agree)
     
      stipple_raster <- raster(resolution=c(0.25, 0.25), ext=extent(agr_data[[p]]))
      
      #Resample so points larger
      fut_mod_agr <- resample(agr_data[[p]], stipple_raster)
      
      
      #Find pixels where model's don't agree
      ind    <- which(values(fut_mod_agr) == 0)
      coords <- coordinates(fut_mod_agr)
      
      #Add stippling
      points(coords[ind,1], coords[ind,2], pch=20, lwd=lwd, cex=cex)
      
      # #Panel number
      # if (p ==1 ){
      #   mtext(side=3, line=0, adj=0, text="c", font=2, xpd=NA, cex=panel_cex) #panel number
      #   
      # } else {
      #   mtext(side=3, line=0, adj=0, text="d", font=2, xpd=NA, cex=panel_cex) #panel number
      #   
      # }
      
      if (v==1 & p==1) mtext(side=2, line=2, text=exp_labels[p])
      
    }
    
  } #metrics
  
  dev.off ()
} #variables





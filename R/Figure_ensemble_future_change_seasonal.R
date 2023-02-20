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
metrics <- c("rel_intensity_by_month", "frequency")

#List seasons
seasons <- c("summer", "autumn", "winter", "spring")

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

lims_diff <- list(duration      = c(-1000, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 1000),
                  intensity     = c(-10000, -30, -25, -20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20, 25, 30, 10000),
                  rel_intensity_by_month = c(-10000, -8, -6, -4, -2, 2, 4, 6, 8, 10000), #c(-100, -40, -30, -20, -10, 10, 20, 30, 40, 100),
                  frequency     = c(-100, -15, -10, -5, -2.5, 2.5, 5, 10, 15, 100))


unit <- c("% points", "% points") # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7



#Loop through metrics
for (m in 1:length(metrics)) {
  
  #Loop through experiments
  for (e in 1:length(exp)) {
    
    
    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_Mean_seasonal_changes_in_", metrics[m], "_",
               percentile, "_", scale, "_", exp[e], ".png"),
        height=8.5, width=8.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.3, 0.4, 0.1))
    
    layout(matrix(c(1:4, 5, 6:9, 5, 10:13, 5), nrow=5), heights=c(1, 1, 1, 1, 0.3,
                                                                  1, 1, 1, 1,
                                                                  1, 1, 1, 1))
    #par(mfcol=c(3, 3))
    par(bty="n")
    
    
    #Loop through variables
    for (v in 1:length(vars)) {
      
      #Loop through seasons
      for (s in 1:length(seasons)) {
        
        
        data_files_hist <- list.files(paste0(path, "/Seasonal_drought_metrics/scale_", 
                                             scale, "/historical/", vars[v]),
                                      pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                      full.names=TRUE, recursive=TRUE)
        
        
        #Read datasets
        data_files_rcp <- list.files(paste0(path, "/Seasonal_drought_metrics/scale_", 
                                            scale, "/", exp[e], "/", vars[v]),
                                     pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                     full.names=TRUE, recursive=TRUE)
        

        #Should find 16 files in each case, check
        if (any(c(length(data_files_hist), length(data_files_rcp)) != 16)) {
          stop("Incorrect number of files found")
        }



        data_hist <- brick(lapply(data_files_hist, raster))
        data_rcp  <- brick(lapply(data_files_rcp, raster))
        
        
        #Need to convert frequency to "percentage of time under drought"
        #Divide by the number of years (in this case 36 as using the mean of two 36-yr periods,
        #i.e. 1970-2005 and 2064-2099)
        #Note something a bit funny going on. 18 drought months during hte historical period
        #but would expect ~16 (15% of 108 months (i.e. 36*3, the number of years multiplied by 3 months in each season))
        if (metrics[m] == "frequency") {
          
          data_hist <- data_hist / (36*3) * 100
          data_rcp  <- data_rcp / (36*3) * 100
          
        }
        
        
        
        #Calculate historical model mean
        ens_median_hist <- calc(data_hist, median)
        
        
        #Calculate future change
        
        future_diff_rcp <- data_rcp - data_hist

        #Ensemble mean of future change
        ens_median_rcp <- calc(future_diff_rcp, median)

        
        #Calculate model agreement
        mod_agr_rcp <- calc(future_diff_rcp, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                                                               weights=rep(1, nlayers(data_rcp))))
      
          
        ################
        ### Plotting ###
        ################
        
        #Set pixels with no model agreement to some crazy values
        #(will plot these in grey)
        mask_value <- -1000000000
        
        ens_median_rcp[mod_agr_rcp == 0] <- mask_value
        
       
      
        ### Future change ###
  
        #Plot
        len <- length(lims_diff[[metrics[m]]])
        
        image(ens_median_rcp, breaks=c(mask_value, mask_value+1, lims_diff[[metrics[m]]][2:len]), 
                                       col=c("grey80", cols_diff(len-1)),
              axes=FALSE, ann=FALSE, asp=1)
        
        
        #Australia outline
        map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
        
        
        #Variable label
        #if(p==1) mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
        
        
        #Calculate and print land area where model changes are robust
        area_agree <- area(mod_agr_rcp)
        area_agree[mod_agr_rcp != 1 | is.na(mod_agr_rcp)] <- NA
        
        area_total <- mask(area(mod_agr_rcp), mod_agr_rcp)
        
        land_area <- (sum(values(area_agree), na.rm=TRUE)  / 
                        sum(values(area_total), na.rm=TRUE)) *100
        
        
        mtext(side=3, line=-3, text=paste0(sprintf("%.1f", land_area), "%"), adj=0.1, cex=0.65)
        
        
        ### Stippling (where models don't agree) ###
        
        # stipple_raster <- raster(resolution=c(0.25, 0.25), ext=extent(mod_agr_rcp))
        # 
        # #Resample so points larger
        # fut_mod_agr <- resample(mod_agr_rcp, stipple_raster)
        # 
        # 
        # #Find pixels where model's don't agree
        # ind    <- which(values(fut_mod_agr) == 0)
        # coords <- coordinates(fut_mod_agr)
        # 
        # #Add stippling
        # points(coords[ind,1], coords[ind,2], pch=20, lwd=lwd, cex=cex)
        # 

        ### Variable label ###
        if (s==1) mtext(side=3, line=1, font=2, text=var_labels[v], xpd=NA)
        
        
        ### Season label ###
        if (v==1) mtext(side=2, line=1, text=seasons[s])
        

      
        ### Legend ###
        
        if (v==1 & s==4) {
          
          #Empty plot
          plot(1, type="n", bty="n", yaxt="n", xaxt="n")
          
          #Legend
          len <- length(lims_diff[[metrics[m]]])-1
          add_raster_legend2(cols=cols_diff(len), limits=lims_diff[[metrics[m]]][2:len],
                             main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                             title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
        }
        
        
  
        
      } #seasons
    } #variables
  
    dev.off ()
    
  } #experiments
  
} #metrics







library(raster)
library(RColorBrewer)
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
metrics <- c("frequency", "rel_intensity_by_month")

#List seasons
seasons <- c("summer", "autumn", "winter", "spring")


outdir <- paste0(path, "/Figures_CABLE_AWRA")
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
  

  #Loop through variables
  for (v in 1:length(vars)) {
    
  
  
    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_Mean_seasonal_changes_CABLE_VS_AWRA_in_", metrics[m], "_",
               percentile, "_", scale, "_", vars[v], ".png"),
        height=10.5, width=8.3, units="in", res=400)
    
    
    par(mai=c(0, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.3, 0.4, 0.1))
    
    layout(matrix(c(1:12, 13, 13, 13), byrow=TRUE, nrow=5), heights=c(1, 1, 1, 1, 0.3,
                                                                  1, 1, 1, 1, 0.3, 
                                                                  1, 1, 1, 1, 0.3))
    
    
    #par(mfcol=c(3, 3))
    par(bty="n")
    
    
       
      #Loop through seasons
      for (s in 1:length(seasons)) {
        
        
        ############
        ### AWRA ###
        ############
        
        data_files_hist_awra <- list.files(paste0(path, "/Seasonal_drought_metrics/scale_", 
                                             scale, "/historical/", vars[v], "/MRNBC/"),
                                      pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                      full.names=TRUE, recursive=TRUE)
        
        
        #Read datasets
        data_files_rcp_awra <- list.files(paste0(path, "/Seasonal_drought_metrics/scale_", 
                                            scale, "/rcp45/", vars[v], "/MRNBC/"),
                                     pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                     full.names=TRUE, recursive=TRUE)
        

        #Should find 16 files in each case, check
        if (any(c(length(data_files_hist_awra), length(data_files_rcp_awra)) != 4)) {
          stop("Incorrect number of files found")
        }

        
        #################
        ### CABLE CO2 ###
        #################
          
        data_files_hist_CO2 <- list.files(paste0(path, "/Seasonal_drought_metrics_CABLE/scale_", 
                                             scale, "/historical/", vars[v], "/CO2/"),
                                      pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                      full.names=TRUE, recursive=TRUE)
        
        
        data_files_rcp_CO2 <- list.files(paste0(path, "/Seasonal_drought_metrics_CABLE/scale_", 
                                                 scale, "/rcp45/", vars[v], "/CO2/"),
                                          pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                         full.names=TRUE, recursive=TRUE)
        
        
        #Should find 16 files in each case, check
        if (any(c(length(data_files_hist_CO2), length(data_files_rcp_CO2)) != 4)) {
          stop("Incorrect number of files found")
        }
        
        
        #################
        ### CABLE noCO2 ###
        #################
        
        data_files_hist_noCO2 <- list.files(paste0(path, "/Seasonal_drought_metrics_CABLE/scale_", 
                                                 scale, "/historical/", vars[v], "/noCO2/"),
                                          pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                          full.names=TRUE, recursive=TRUE)
        
        
        data_files_rcp_noCO2 <- list.files(paste0(path, "/Seasonal_drought_metrics_CABLE/scale_", 
                                                scale, "/rcp45/", vars[v], "/noCO2/"),
                                         pattern=glob2rx(paste0("Mean_", metrics[m], "*", seasons[s], ".nc")),
                                         full.names=TRUE, recursive=TRUE)
        
        
        #Should find 16 files in each case, check
        if (any(c(length(data_files_hist_noCO2), length(data_files_rcp_noCO2)) != 4)) {
          stop("Incorrect number of files found")
        }
        
        
        
        
        #################
        ### Load data ###
        #################
        
        
        data_hist_awra  <- brick(lapply(data_files_hist_awra, raster))
        data_rcp_awra   <- brick(lapply(data_files_rcp_awra, raster))
        
        data_hist_CO2   <- brick(lapply(data_files_hist_CO2, raster))
        data_rcp_CO2    <- brick(lapply(data_files_rcp_CO2, raster))
        
        data_hist_noCO2 <- brick(lapply(data_files_hist_noCO2, raster))
        data_rcp_noCO2  <- brick(lapply(data_files_rcp_noCO2, raster))
        
        
        
        #Need to convert frequency to "percentage of time under drought"
        #Divide by the number of years (in this case 36 as using the mean of two 36-yr periods,
        #i.e. 1970-2005 and 2064-2099)
        #Note something a bit funny going on. 18 drought months during hte historical period
        #but would expect ~16 (15% of 108 months (i.e. 36*3, the number of years multiplied by 3 months in each season))
        if (metrics[m] == "frequency") {
          
          data_hist_awra  <- data_hist_awra / (36*3) * 100
          data_rcp_awra   <- data_rcp_awra / (36*3) * 100
          
          data_hist_CO2   <- data_hist_CO2 / (36*3) * 100
          data_rcp_CO2    <- data_rcp_CO2 / (36*3) * 100
          
          data_hist_noCO2 <- data_hist_noCO2 / (36*3) * 100
          data_rcp_noCO2  <- data_rcp_noCO2 / (36*3) * 100

        }
        
        
        
        #Calculate future change
        future_diff_rcp_awra  <- data_rcp_awra - data_hist_awra
        future_diff_rcp_CO2   <- data_rcp_CO2 - data_hist_CO2
        future_diff_rcp_noCO2 <- data_rcp_noCO2 - data_hist_noCO2
        
        
        #Ensemble mean of future change
        ens_median_awra  <- calc(future_diff_rcp_awra, median)
        ens_median_CO2   <- calc(future_diff_rcp_CO2, median)
        ens_median_noCO2 <- calc(future_diff_rcp_noCO2, median)
        
        
      
          
        ################
        ### Plotting ###
        ################
        
         
        ### Model differences ###
        
        plot_data <- list(awra_CO2   = ens_median_awra - ens_median_CO2, #(ens_median_awra - ens_median_CO2) / ens_median_CO2 *100,
                          awra_noCO2 = ens_median_awra - ens_median_noCO2,
                          CO2_noCO2  = ens_median_noCO2 - ens_median_CO2)
        
        
        labs <- c("AWRA - CO2", "AWRA - noCO2", "noCO2 - CO2")
        
        cols <- colorRampPalette(rev(c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", 
                                       "#c7eae5", "#80cdc1", "#35978f", "#01665e")))
                                                
        #Loop through experiments
        for (p in 1:length(plot_data)) {
          
          
          #Plot
          len <- length(lims_diff[[metrics[m]]])
          image(plot_data[[p]], breaks=lims_diff[[metrics[m]]], 
                col=c(cols(len-1)),
                axes=FALSE, ann=FALSE, asp=1)
          
          
          #Australia outline
          map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
          
          
          #Model label
          if (s==1) mtext(side=3, line=0, font=2, text=labs[p], xpd=NA)
          
          
          if (p==1) mtext(side=2, line=1, text=seasons[s])
          
          
        }
        
        
        if (s == length(seasons)) {
          
          ### Legend ###
          
          #Empty plot
          plot(1, type="n", bty="n", yaxt="n", xaxt="n")
          
          #Legend
          len1 <- length(lims_diff[[metrics[m]]])-1
          add_raster_legend2(cols=cols(len1), limits=lims_diff[[metrics[m]]][2:len1],
                             main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                             title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
          
        }
         
        
  
        
      } #seasons
    
      dev.off ()
    } #variables
    
} #metrics







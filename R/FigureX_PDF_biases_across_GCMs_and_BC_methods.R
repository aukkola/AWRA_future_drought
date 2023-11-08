library(raster)
library(RColorBrewer)
library(maptools)
library(maps)
library(parallel)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"


source(paste0(path, "/scripts/R/functions/trend_per_pixel.R"))
source(paste0(path, "/scripts/R/functions/pdf_and_skill.R"))

#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3


#Variables
vars <- c("pr", "qtot", "sm")#, "mrro") #list.files(paste0(dr_path, exp[1]))

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")


#Experiments
exp <- c("rcp45")#, "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Output directory
outdir <- paste0(path, "/Figures")
dir.create(outdir)



#Set up parallel processing for trend calculation
cl <- makeCluster(getOption('cl.cores', 16))

clusterExport(cl, c('trend_per_pixel_MannKen', 'trend_per_pixel'))
clusterEvalQ(cl, library(raster))



#####################
### Plot settings ###
#####################


cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")


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



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  ### Set up figure ###
  png(paste0(outdir, "/FigureX", "_PDFs_of_gmcs_and_bc_methods_in_", metrics[m], "_",
             percentile, "_", scale, "_", exp, ".png"),
      height=9.5, width=6.3, units="in", res=400)
  
  
  par(mai=c(0, 0.1, 0.2, 0.1))
  par(omi=c(0.1, 0.4, 0.4, 0.3))
  
  par(mfrow=c(length(vars), 2))
  
  # layout(matrix(c(1:4, 5, 6:9, 5, 10:13, 5), nrow=5), heights=c(1 ,0.3, 1, 1, 0.3,
  #                                                               1 ,0.3, 1, 1, 0.3,
  #                                                               1 ,0.3, 1, 1, 0.3))
  # #par(mfcol=c(3, 3))
  #par(bty="n")
  
  

  #Loop through variables
  for (v in 1:length(vars)) {
  
     
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
    data <- lapply(data_files, brick, varname=metrics[m])
    
    trend <- parLapply(cl=cl, data, function(x) calc(x, trend_per_pixel)) #_MannKen))
   
    
    # #Calculate mean trend across all simulations
    # mean_trend <- mean(brick(lapply(trend, function(x) x[[1]])))
    # 
    # 
    # #Then calculate percentage difference from mean trend
    # rel_trend <- lapply(trend, function(x) (x - mean_trend)/x*100)
    # 
    
    
    ################
    ### Plotting ###
    ################
    
    #Function to calculate probability densities
    hist_xy <- function(data, n=50) {

      hist <- hist(data, breaks=n, plot=FALSE)
   #   offset <- hist$breaks[2] - hist$breaks[1]
    #  density <- hist$density*offset #normalise density to sum up to 1
      return(list(x=hist$breaks, y=hist$density))
    }
    
    
    #Calculate area of each pixel (make sure to mask as area function returns area for all pixels, including NA)
    area <- mask(area(trend[[1]][[1]]), trend[[1]][[1]])
    
    #Determine range of trend values to set histogram bins
    range <- range(unlist(trend, function(x) range(values(x), na.rm=TRUE)))
    
    #Histogram breaks
    breaks <- seq(range[1], range[2], length.out=100)
    
    
    
    #######################
    ### Plot PDF by GCM ###
    #######################
    
    gcm_ind <- lapply(gcms, function(x) which(grepl(x, data_files)))
    
    #Get trend values 
    hist_data_gcm <- lapply(1:length(gcms), function(x) lapply(gcm_ind[[x]], function(y) values(trend[[y]][[1]])))
    
    hist_gcm <- lapply(hist_data_gcm, function(x) freqs(unlist(x), breaks, area))
    
    #Get plot ranges
    plot_y_range <- range(unlist(lapply(hist_gcm, function(x) x)), na.rm=TRUE)
    
    #x-range has some crazy values from the percentage calculation, ignore these
    #by taking x-range as the 10/90th percentile
    plot_x_range <- quantile(unlist(lapply(hist_gcm, function(x) x$x)), 
                             probs=c(0.1, 0.9), na.rm=TRUE)
    
    
    #Plotting
    plot(plot_x_range, plot_y_range, type="n", ylab="n", xlab="n")
    
    
    #Add model lines
    for (gcm in 1:length(hist_gcm)) {
      
      #Get x data and center it (hist function returns min and max of each bin,
      #want to use the middle value)
      
      x <- hist_gcm[[gcm]]$x
      
      offset <- x[2] - x[1]
      
      x_centred <- x[2:length(x)] - offset
      
      #Plot line
      lines(x_centred, hist_gcm[[gcm]]$y, col=cols[gcm])
    }
    
    if(v==1) {
      legend("topright", legend=gcm_labels, col=cols, bty="n", lty=1)
      mtext(side=2, line=2, var_labels[v])
      mtext(side=3, line=1, "GCM")
      
    }
    
    #############################
    ### Plot PDF by BC method ###
    #############################
    
    bc_ind <- lapply(bc_methods, function(x) which(grepl(x, data_files)))
    
    #Get trend values 
    hist_data_bc <- lapply(1:length(bc_methods), function(x) lapply(bc_ind[[x]], function(y) values(trend[[y]][[1]])))
    
    hist_bc <- lapply(hist_data_bc, function(x) hist_xy(unlist(x)))
    
    #Get plot ranges
    plot_y_range <- range(unlist(lapply(hist_bc, function(x) x$y)), na.rm=TRUE)
    
    #x-range has some crazy values from the percentage calculation, ignore these
    #by taking x-range as the 10/90th percentile
    plot_x_range <- quantile(unlist(lapply(hist_bc, function(x) x$x)), 
                             probs=c(0.1, 0.9), na.rm=TRUE)
    
    
    #Plotting
    plot(plot_x_range, plot_y_range, type="n", ylab="n", xlab="n")
    
    
    #Add model lines
    for (bc in 1:length(hist_bc)) {
      
      #Get x data and center it (hist function returns min and max of each bin,
      #want to use the middle value)
      
      x <- hist_bc[[bc]]$x
      
      offset <- x[2] - x[1]
      
      x_centred <- x[2:length(x)] - offset
      
      #Plot line
      lines(x_centred, hist_bc[[bc]]$y, col=cols[bc])
    }
    
    if(v==1) {
      legend("topright", legend=bc_methods, col=cols, bty="n", lty=1)
      mtext(side=3, line=2, "BC method")
    }
    
    
   
  } #variables
  
  dev.off ()
} #metrics


stopCluster(cl)



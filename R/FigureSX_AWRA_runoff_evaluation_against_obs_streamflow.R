library(raster)
library(maptools)
library(parallel)
library(maps)

path <- "/g/data/w97/amu561/Steven_CABLE_runs/"

scratch_path <- "/scratch/w97/amu561/monthly_sums"

source(paste0(path,"/scripts/R/functions/drought_metrics.R"))
source(paste0(path, "/scripts/R/functions/pdf_and_skill.R"))

#Load streamflow data
streamflow_file <- paste0(path, "/../../mg5624/RF_project/Streamflow/processed/monthly", 
                          "_streamflow_data_processed_1970-2020_missing_cutoff_0.05.csv")

streamflow <- read.csv(streamflow_file)

#Load basin shapefiles
shapefile <- paste0(path, "/../../mg5624/RF_project/Streamflow/02_location_boundary_area/", 
                    "shp/CAMELS_AUS_v2_Boundaries_adopted.shp")

basins <- readShapePoly(shapefile)


#Set drought metric settings
scale <- 3
perc <- 0.15


metrics <- c("duration", "intensity")

cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

metric_labels <- c("Duration (months)", "Relative intensity (%)")



################################
### Observed drought metrics ###
################################

#Calculate observed drought duration and intensity

obs_drought_metrics <- lapply(1:(ncol(streamflow)-1), function(x) drought_metrics(streamflow[,x+1], scale=scale, perc=perc))



############################
### AWRA drought metrics ###
############################


#Read datasets
data_files_rcp45 <- list.files(paste0(scratch_path, "/"),
                               pattern=paste0("_rcp45_qtot.nc"),
                               full.names=TRUE)



#Should find 16 files in each case, check
if (length(data_files_rcp45) != 16) {
  stop("Incorrect number of files found")
}



data_rcp45 <- lapply(data_files_rcp45, brick)

#Expected number of months
nmonths <- length(1960:2099)*12

if(any(sapply(data_rcp45, nlayers) != nmonths)) {
  stop("incorrect number of months") }

#Get time stamps
dates <- getZ(data_rcp45[[1]])

#Get start and end dates matching Jan 1970, Dec 2020
start <- which(grepl("1970", dates))[1]
end   <- tail(which(grepl("2020", dates)), 1)

#Get correct time period
data_rcp45 <- lapply(data_rcp45, function(x) x[[start:end]])



#Get IDs of basins that have data during 1970-2020 (first column is the date, remove)

basin_IDs <- colnames(streamflow)[2:ncol(streamflow)]

#Remove extra "X" character that R inserts when reading the csv
basin_IDs <- gsub("X", "", basin_IDs)

#List all available basins
all_basins <- basins$CatchID


#intersect and find indices for available basins
basin_indices <- sapply(basin_IDs, function(x) which(all_basins == x))

#Check that found all
if(length(basin_indices) != length(basin_IDs)) stop("Some basins not found")


#Save output when created, takes a long time to run

basin_outputfile <- paste0(path, "/AWRA_obs_streamflow_comparison/AWRA_basin_runoff.rds")

if (!file.exists(basin_outputfile)) {
  
  #Use parallel to calculate 3-month running means
  cl <- makeCluster(getOption('cl.cores', 16))
  
  clusterExport(cl, c('basins', 'basin_indices', 'extract'))
  
  awra_basin_data <- parLapply(cl=cl, data_rcp45, function(x) extract(x, basins[basin_indices,], 
                                                                      fun=mean, weights=TRUE))
  
  stopCluster(cl)
  
  dir.create(paste0(path, "/AWRA_obs_streamflow_comparison/"))
  
  saveRDS(awra_basin_data, basin_outputfile)
  
} else {
  
  awra_basin_data <- readRDS(basin_outputfile)
  
}




### Split by GCM/BC-DS method as per Figure 5 in main paper ###

#Get BC methods
bc_methods <- list.files(paste0(path, "/drought_metrics/", 
                                scale, "-month/"))

#Get GCMs
gcms <- list.files(paste0(path, "/drought_metrics/", 
                          scale, "-month/CCAM/"))


#Prettier labels
gcm_labels <- c(CNRM  = "CNRM-CM5",
                CSIRO = "ACCESS1-0",
                MIROC = "MIROC5",         
                NOAA  = "GFDL-ESM2M")


### GCM ###

#Find files for each GCM
gcm_ind <- lapply(gcms, function(x) which(grepl(x, data_files_rcp45)))

### BC ###

#Find files for each BC-DS method
bc_ind <- lapply(bc_methods, function(x) which(grepl(x, data_files_rcp45)))



###############
### Metrics ###
###############


#Function to calculate frequencies
freqs <- function(x, breaks) {
  h <- hist(x, breaks=breaks, plot=FALSE)
  density <- h$density #this is the same ash$counts / (sum(h$counts) * diff(h$breaks)) 
  #and should be called probability density because it adjusts for bin width. A density 
  #represents the probability per unit on the x-axis, and its integral (area under the curve) 
  #over the entire range should equal 1. To calculate density properly for a histogram, 
  #you need to account for the width of each bin because the probability density scales 
  #with both bin counts and bin widths.
  
  #density <- h$counts / sum(h$counts) #this gives the relative frequency
  return(density)
}


### Set up figure ###
png(paste0(path, "/Figures/FigureSX_historical_PDFs_in_runoff_for_observed_catchments.png"),
    height=6.5, width=6.3, units="in", res=400)


par(mai=c(0.6, 0.1, 0.2, 0.1))
par(omi=c(0.5, 0.5, 0.4, 0.3))

par(mfrow=c(length(metrics), 2))



for (m in 1:length(metrics)) {
  
  #Get observed metric
  obs_metric <- unlist(sapply(obs_drought_metrics, function(x) x[metrics[m]]))
  
  #Calculate model metrics
  awra_metric <- lapply(awra_basin_data, function(x) unlist(apply(x, function(y) 
                        drought_metrics(y, scale, perc)[metrics[m]], MARGIN=1)))
  
  
  
  
  ### Set up histogram settings ###
  
  #Number of bins
  nbins <- 15
  
  data_range <- range(obs_metric, unlist(awra_metric), na.rm=TRUE)
  
  breaks <- seq(data_range[1], data_range[2], length.out=nbins)
  

  
  
  ### Calculate densities ###
  
  gcm_freqs <- lapply(gcm_ind, function(x) freqs(unlist(awra_metric[x]), breaks))

  bc_freqs <- lapply(bc_ind, function(x) freqs(unlist(awra_metric[x]), breaks))
  
  obs_freqs <- freqs(obs_metric, breaks)
  
  
  #Plot ranges
  plot_y_range <- range(unlist(gcm_freqs), bc_freqs, obs_freqs, na.rm=TRUE)
  plot_x_range <- range(breaks)
  
  
  ### GCMs ###
  
  #Plot empty plot
  plot(plot_x_range, plot_y_range, type="n", ylab="", xlab="")
  
  
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
  
  if(m==1) {
    mtext(side=3, line=1, cex=1.2, "GCM")
    
  }
  
  #Metric label
  #mtext(side=2, line=5, cex=1.2, metrics[m])
  
  #y-label
  mtext(side=2, line=2.5, "Probability density (-)")
  
  #x-label
  mtext(side=1, line=2.5, metric_labels[m])

  
  ### Add legend ###
  
  
  #Add legend with skill scores
  perkins_gcm <- sapply(gcm_freqs, function(x) perkins_skill_score_normalised(x, obs_freqs))
  mbe_gcm     <- sapply(gcm_ind, function(x) mean(unlist(awra_metric[x]), na.rm=TRUE) - mean(obs_metric, na.rm=TRUE))
  
  
  #Need to place legend in a different spot depending on variable/metric
  position <- "topright"
  
  if (metrics[m] == "intensity" ) position <- "topleft"
  
  legend(position, legend=paste0(gcm_labels, " (s=", round(perkins_gcm, digits=2), ", MBE=", 
                                 round(mbe_gcm, digits=2), ")"), col=cols, bty="n", lty=1, cex=0.7)
  
  
  
  #####################
  ### BC-DS methods ###
  #####################
  
  
  #Plot empty plot
  plot(plot_x_range, plot_y_range, type="n", ylab="", xlab="", yaxt="n")
  
  
  #Add obs line
  lines(x_centred, obs_freqs, col="black", lwd=1)
  
  #Add model lines
  for (bc in 1:length(bc_freqs)) {
    
    #Plot line
    lines(x_centred, bc_freqs[[bc]], col=cols[bc])
  }
  
  #main title
  if(m==1) {
    mtext(side=3, line=1, cex=1.2, "DS-BC method")
    
  }
  
  #x-label
  mtext(side=1, line=2.5, metric_labels[m])
  
  
  
  ### Add legend ###
  
  #Add legend with skill scores
  perkins_bc <- sapply(bc_freqs, function(x) perkins_skill_score_normalised(x, obs_freqs))
  mbe_bc     <- sapply(bc_ind, function(x) mean(unlist(awra_metric[x]), na.rm=TRUE) - mean(obs_metric, na.rm=TRUE))
  
  #Legend
  legend(position, legend=paste0(bc_methods, " (s=", round(perkins_bc, digits=2), ", MBE=", 
                                 round(mbe_bc, digits=2), ")"), col=cols, bty="n", lty=1, cex=0.7)
  
  
} #metrics



dev.off()





##################################
### Map of included catchments ###
##################################


png(paste0(path, "/Figures/FigureSX_CAMELS_catchments_for_1970-2020.png"),
    height=4, width=4, units="in", res=400)


par(mai=c(0, 0, 0, 0))
par(omi=c(0.1, 0.1, 0.1, 0.1))


#Australia outline
map(region="Australia", lwd=0.7, fill=TRUE, col="grey95", border="grey50")


#Plot basin boundaries
plot(basins[basin_indices,], col="#66c2a5", add=TRUE, lwd=0.2)

#text(120, -40, "n = 216", cex=1)


dev.off()






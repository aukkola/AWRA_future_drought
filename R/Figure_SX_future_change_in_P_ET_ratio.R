library(raster)
library(RColorBrewer)
library(maptools)
library(maps)
library(zoo)
library(parallel)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"

scratch_path <- "/scratch/w97/amu561/monthly_sums"

#Source functions
source(paste0(path,"/scripts/R/functions/perc_agree_on_sign_weighted_AR4_method.R"))
source(paste0(path,"/scripts/R/functions/wtd_stdev.R"))
source(paste0(path,"/scripts/R/functions/add_raster_legend.R"))


#Load world shapefile
world <- readShapePoly(paste0(path, "../World_shapefile/World"))

scale <- 3


#Variables
vars <- c("pr", "etot")

var_labels <- c("Precipitation", "Evapotranspiration") #labels for plotting


#Experiments
exp <- c("rcp45", "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Raster 1x1 deg for resampling to common resolution
#extent_raster  <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)


#####################
### Set up figure ###
#####################

outdir <- paste0(path, "/Figures")
dir.create(outdir)


col_opts <- rev(brewer.pal(11, 'RdYlBu'))
#col_opts[6] <- "grey80"
cols_diff <- colorRampPalette(col_opts) 

#Limits

lims_diff <- c(-1000, -0.05, -0.04, -0.03, -0.02, 0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 1000)


unit <- "(-)" # (#expression("mm"~"month"^"-1"), expression("no. events 10 yrs"^"-1"))

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7


### Set up figure ###
png(paste0(outdir, "/FigureSX", "_Ensemble_median_change_in_ET_pr_ratio", ".png"),
    height=4, width=7.5, units="in", res=400)


par(mai=c(0, 0.0, 0.1, 0.0))
par(omi=c(0.1, 0.5, 0.4, 0.1))

#par(mfcol=c(3, 3))
par(bty="n")


#layout(matrix(c(1:4, 5, 6:9, 5), nrow=5), heights=c(1 , 1, 0.3))

layout(matrix(c(1, 2, 3, 3), nrow=2, byrow=TRUE), heights=c(1 , 0.3))

#Initialise
hist_all_data  <- list()
rcp45_all_data <- list()
rcp85_all_data <- list()


#Loop through variables
for (v in 1:length(vars)) {
  
  #Read datasets
  data_files_rcp45 <- list.files(paste0(scratch_path, "/"),
                                 pattern=paste0("_rcp45_", vars[v], ".nc"),
                                 full.names=TRUE)
  
  data_files_rcp85 <- list.files(paste0(scratch_path, "/"),
                                 pattern=paste0("_rcp85_", vars[v], ".nc"),
                                 full.names=TRUE)
  
  #Should find 16 files in each case, check
  if (any(c(length(data_files_rcp45),
            length(data_files_rcp85)) != 16)) {
    stop("Incorrect number of files found")
  }
  
  
  
  data_rcp45 <- lapply(data_files_rcp45, brick)
  data_rcp85 <- lapply(data_files_rcp85, brick)
  
  #Expected number of months
  nmonths <- length(1960:2099)*12
  
  if(any(c(sapply(data_rcp45, nlayers), sapply(data_rcp85, nlayers)) != nmonths)) {
    stop("incorrect number of months") }
  
  #Calculate the difference over:
  #1970-2005
  #2064-2099
  
  #Can't use getZ to obtain dates for some reason. Check that both rcp45 and rcp85
  #datasets have the expected number of time slices and extract the time periods 
  #manually
  
  hist_start <- 121 #Jan 1970
  hist_end   <- 552 #Dec 2005
  
  fut_start <- 1249 #Jan 2064
  fut_end   <- 1680 #Dec 2099
  
  
  #Use parallel to calculate 3-month running means
  cl <- makeCluster(getOption('cl.cores', 16))
  
  clusterExport(cl, c('rollmean', 'scale', 'calc', 'hist_start', 'hist_end',
                      'fut_start', 'fut_end', 'mean', 'sd'))
  #clusterEvalQ(cl, library(raster))
  
  
  
  # # Rolling means
  # data_rcp45_rollmean <- parLapply(cl=cl, data_rcp45, function(x) calc(x, function(y) 
  #   rollmean(y, scale, na.pad=TRUE)))
  # data_rcp85_rollmean <- parLapply(cl=cl, data_rcp85, function(x) calc(x, function(y) 
  #   rollmean(y, scale, na.pad=TRUE)))
  # 
  
  ### Future and historical mean ###
  
  #Calculate historical and future mean pr/ET (only need to calculate one historical mean,
  #the same across both RCPs)
  hist_mean <- brick(parLapply(cl=cl, data_rcp45, function(x) mean(x[[hist_start:hist_end]])))
  
  rcp45_fut_mean <- brick(parLapply(cl=cl, data_rcp45, function(x) mean(x[[fut_start:fut_end]])))
  rcp85_fut_mean <- brick(parLapply(cl=cl, data_rcp85, function(x) mean(x[[fut_start:fut_end]])))
  
  stopCluster(cl)
  
  
  #Collate
  hist_all_data[[vars[v]]]  <- hist_mean
  rcp45_all_data[[vars[v]]] <- rcp45_fut_mean
  rcp85_all_data[[vars[v]]] <- rcp85_fut_mean
  
  
  
} #vars
 

################
### Plotting ###
################

##########################
### ET/P ratio changes ###
##########################

#Calculate future change

hist_ratio <- hist_all_data[["etot"]] / hist_all_data[["pr"]]

rcp45_ratio <- rcp45_all_data[["etot"]] / rcp45_all_data[["pr"]]
rcp85_ratio <- rcp85_all_data[["etot"]] / rcp85_all_data[["pr"]]

fut_diff_rcp45 <- rcp45_ratio - hist_ratio
fut_diff_rcp85 <- rcp85_ratio - hist_ratio

#Ensemble mean of future change
ens_median_rcp45 <- calc(fut_diff_rcp45, median)
ens_median_rcp85 <- calc(fut_diff_rcp85, median)


#Calculate model agreement
mod_agr_rcp45 <- calc(fut_diff_rcp45, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                                                          weights=rep(1, nlayers(fut_diff_rcp45))))
mod_agr_rcp85 <- calc(fut_diff_rcp45, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                                                          weights=rep(1, nlayers(fut_diff_rcp45))))

#List and then loops
plot_data <- list(rcp45 = ens_median_rcp45,
                  rcp85 = ens_median_rcp85)

agr_data <- list(rcp45 = mod_agr_rcp45,
                 rcp85 = mod_agr_rcp85)



#Loop through experiments
for (p in 1:length(plot_data)) {
  
  #Set pixels with no model agreement to some crazy values
  #(will plot these in grey)
  # mask_value <- -1000000000
  
  #  plot_data[[p]][agr_data[[p]] == 0] <- mask_value
  
  
  #Plot
  len <- length(lims_diff)
  image(plot_data[[p]], breaks=lims_diff, 
        col=c(cols_diff(len-1)),
        axes=FALSE, ann=FALSE, asp=1)
  
  
  #Australia outline
  map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
  
  
  #Variable label
  mtext(side=3, line=1, font=2, text=exp_labels[p], xpd=NA)
  
  
  # #Calculate and print land area where model changes are robust
  # area_agree <- area(agr_data[[p]])
  # area_agree[agr_data[[p]] != 1 | is.na(agr_data[[p]])] <- NA
  # 
  # area_total <- mask(area(agr_data[[p]]), agr_data[[p]])
  # 
  # land_area <- (sum(values(area_agree), na.rm=TRUE)  / 
  #                 sum(values(area_total), na.rm=TRUE)) *100
  # 
  # 
  # mtext(side=3, line=-3, text=paste0(sprintf("%.1f", land_area), "%"), adj=0.1, cex=0.65)
  
  
}

#Legend

  
#Empty plot
plot(1, type="n", bty="n", yaxt="n", xaxt="n")

#Legend
len <- length(lims_diff)-1
add_raster_legend2(cols=cols_diff(len), limits=lims_diff[2:len],
                   main_title=paste0("Change in ET/P ratio (-)"), 
                   plot_loc=c(0.3,0.7,0.35, 0.5), 
                   title.cex=1.1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE,
                   title_fac=0.2)




dev.off()
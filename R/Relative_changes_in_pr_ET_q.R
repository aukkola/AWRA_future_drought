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


scale <- 3


#Variables
vars <- c("pr", "etot", "qtot")


#Experiments
exp <- c("rcp45", "rcp85")


#####################
### Set up figure ###
#####################

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


#Initialise
hist_all_data  <- list()
rcp45_all_data <- list()


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

  #Expected number of months
  nmonths <- length(1960:2099)*12
  
  if(any(c(sapply(data_rcp45, nlayers)) != nmonths)) {
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

  stopCluster(cl)
  
  
  #Collate
  hist_all_data[[vars[v]]]  <- hist_mean
  rcp45_all_data[[vars[v]]] <- rcp45_fut_mean
  
  
} #vars


pr_ratio <- (rcp45_all_data[["pr"]] - hist_all_data[["pr"]]) / hist_all_data[["pr"]] *100
et_ratio <- (rcp45_all_data[["etot"]] - hist_all_data[["etot"]]) / hist_all_data[["etot"]] *100
q_ratio <-  (rcp45_all_data[["qtot"]] - hist_all_data[["qtot"]]) / hist_all_data[["qtot"]] *100


pr_ratio[agr_data[[1]] < 1] <- NA
et_ratio[agr_data[[1]] < 1] <- NA
q_ratio[agr_data[[1]] < 1] <- NA

pr_ratio[pr_ratio < 0 & q_ratio > 0] = NA
pr_ratio[pr_ratio > 0 & q_ratio < 0] = NA

et_ratio[et_ratio > 0 & q_ratio > 0] = NA
et_ratio[et_ratio < 0 & q_ratio < 0] = NA




par(mfcol=c(2,2))
plot(calc(pr_ratio, median))

plot(calc(et_ratio, median))

plot(calc(q_ratio, median))

pr_median = calc(pr_ratio, median)
et_median = calc(et_ratio, median)
q_median = calc(q_ratio, median)


pr_vals <- values(pr_ratio) #values(pr_median)
et_vals <- values(et_ratio) #values(et_median)
q_vals <- values(q_ratio) #values(q_median)


mean(abs(pr_vals), na.rm=TRUE)
mean(abs(et_vals), na.rm=TRUE)
mean(abs(q_vals), na.rm=TRUE)




ratio_to_use <- et_ratio
q_ratio[is.na(ratio_to_use)] = NA
test= abs(q_ratio) > abs(ratio_to_use)

frac <-vector()
for(n in 1:nlayers(test)) {
  
  vals <- values(test[[n]])
  
  vals <- vals[-which(is.na(vals))]
  
  frac[n] <- round(length(which(vals==1)) / length(vals), 2)
  
}









##################################################
# Drought metric for masking

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
vars <- c("qtot", "pr")#, "mrro") #list.files(paste0(dr_path, exp[1]))



#List metrics
metrics <- c("frequency")


#Experiments
exp <- c("rcp45", "rcp85")

#Level of model agreement for stippling (as fraction of models)
agr_level <- 0.75


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
    
    
    data_hist  <- brick(lapply(data_files_hist, raster))
    data_rcp45 <- brick(lapply(data_files_rcp45, raster))

    
    #Calculate historical model mean
    ens_median_hist <- calc(data_hist, median)
    
    
    #Calculate future change
    
    future_diff_rcp45 <- data_rcp45 - data_hist

    #Ensemble mean of future change
    ens_median_rcp45 <- calc(future_diff_rcp45, median)

    
    #Calculate model agreement
    mod_agr_rcp45 <- calc(future_diff_rcp45, function(x) perc_agree_on_sign_weighted_AR4_method(x, 
                                                                                                weights=rep(1, nlayers(data_rcp45))))

    
    agr_data <- list(rcp45 = mod_agr_rcp45)
    
  }} 
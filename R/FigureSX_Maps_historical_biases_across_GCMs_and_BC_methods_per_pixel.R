library(raster)
library(RColorBrewer)
library(maptools)
library(maps)
library(parallel)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" #"/srv/ccrc/data04/z3509830/CMIP6_drought//"

source(paste0(path, "/scripts/R/functions/pdf_and_skill.R"))
source(paste0(path, "/scripts/R/functions/add_raster_legend.R"))


#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3


#Variables
vars <- c("pr", "qtot", "sm_root")

var_labels <- c("Precipitation", "Runoff", "Soil moisture") #labels for plotting


#List metrics
metrics <- c("duration", "rel_intensity")#, "frequency")

metric_labels <- c("Duration (months)", "Relative intensity (%)")#, "frequency") #labels for plotting


#Experiments
exp <- c("rcp45")#, "rcp85")

exp_labels <- c("RCP-4.5", "RCP-8.5")

#Output directory
outdir <- paste0(path, "/Figures")
dir.create(outdir)



#####################
### Plot settings ###
#####################


cols <- colorRampPalette(rev(c("#b2182b", "#d6604d", "#f4a582", "#fddbc7",
                               "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac")))
                                            

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


plot_range <- list(duration=c(-10000, seq(-1, 1, by=0.2), 10000),
                   rel_intensity=c(-10000, seq(-25,25, by=5), 1000))

#Set up parallel processing for trend calculation
#cl <- makeCluster(getOption('cl.cores', 28))

#clusterEvalQ(cl, library(raster))
beginCluster(28)

#Plot GCM and BC in separate plot
#plot_type <- c("GCM", "BC")


#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  #Progress
  print(paste0("plotting metric ", m, "/", length(metrics)))
  
  # #Plot type (GCM/BC)  
  # for (p in 1:length(plot_type)) {
  #   
  
   
  
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
   
    ### Set up figure ###
    png(paste0(outdir, "/FigureS6-8", "_historical_map_of_biases_in_", metrics[m], "_",
               percentile, "_", scale, "_", vars[v], ".png"),
        height=9.5, width=5.3, units="in", res=400)
    
    
    par(mai=c(0.1, 0.1, 0.1, 0.1))
    par(omi=c(0.1, 0.1, 0.6, 0.1))
    
    layout(matrix(c(1:4, 9, 5:8,9), nrow=5), heights=c(1 ,1, 1, 1, 0.3,
                                                     1 ,1, 1, 1, 0.3))
    #par(mfcol=c(3, 3))
    par(bty="n")
    
    
    print(paste0("variable ", v, "/", length(vars)))
    
    
    ################
    ### Obs data ###
    ################
    
    if (vars[v] == "pr") {
      
      obs_file <- paste0(path, "/drought_metrics_AGCD/", scale, "-month/",
                         "drought_metrics_AGCD_precip_1900_2021_baseline_1970_2005_scale_3.nc")
      
      #Select years 1970-2020 to match model runs and AWRA reference runs
      obs_data <- brick(obs_file, varname=metrics[m])[[841:1452]]
      
    } else {
      
      #Use different variable name for sm. should recreate file but leave for now
      variable <- vars[v]
      if(vars[v] == "sm_root"){
        variable <- "sm"
      } 
      
      obs_file <- paste0(path, "/drought_metrics_AWRA_ref/", scale, "-month/",
                         "drought_metrics_AWRA_ref_", variable, "_scale_", scale, "_1960_2020.nc")
      
      #Select 1970-2020
      obs_data <- brick(obs_file, varname=metrics[m])
      obs_data <- obs_data[[121:nlayers(obs_data)]]
    }
    
    #Calculate obs mean
    obs_mean <- clusterR(obs_data, mean, args=list(na.rm=TRUE))
    
    
    
    ##################
    ### Model data ###
    ##################
    
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
    #Only get data for 1970-2020 to match observations (model data runs 1960-2099)
    data <- lapply(data_files, function(x) brick(x, varname=metrics[m])[[121:732]])
    
    start=Sys.time()
    #Calculate model means
    model_mean <- lapply(data, function(x) clusterR(x, mean, args=list(na.rm=TRUE))) #parLapply(cl, data, mean, na.rm=TRUE)
    end=Sys.time()
    print(paste("time: ", end-start))
    
    
    #Mask rainfall obs with models (obs includes ocean areas)
    if (vars[v]=="pr") {
      
      obs_mean <- mask(crop(obs_mean, model_mean[[1]]),  model_mean[[1]])
    } 
    
    
    ################
    ### Plotting ###
    ################
 
    ########################
    ### Calculate biases ###
    ########################
    
    #Calculate difference to obs
    
    #GCM
    gcm_ind <- lapply(gcms, function(x) which(grepl(x, data_files)))
    
    gcm_biases <- lapply(gcm_ind, function(x) mean(brick(model_mean[x])) - obs_mean)
    
    #BC
    bc_ind <- lapply(bc_methods, function(x) which(grepl(x, data_files)))
    
    bc_biases <- lapply(bc_ind, function(x) mean(brick(model_mean[x])) - obs_mean)
    
    
    
    ###################
    ### Plot biases ###
    ###################
    
    #Collate
    plot_data <- append(gcm_biases, bc_biases)
    
    #Plot labels (GCM or BC name)
    main_labels <- c(gcm_labels, bc_methods)
    
    #breaks
    breaks <- plot_range[[metrics[m]]]
    
    for(p in 1:length(plot_data)) {
      
      #Plot
      image(plot_data[[p]], col=cols(length(breaks)-1), breaks=breaks,
            axes=FALSE, ann=FALSE, asp=1)
      
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
      
      
      #main title
      mtext(side=3, line=0, main_labels[p], cex=0.8)
      
      if(p==1)  mtext(side=3, line=3, "GCM", cex=1.2) #GCM label
      if(p==5)  mtext(side=3, line=3, "DS-BC method", cex=1.2) #GCM label
      
    }
    
    
    #Empty plot
    plot(1, type="n", bty="n", yaxt="n", xaxt="n")
    
    #Legend
    add_raster_legend2(cols=cols(length(breaks)-1), limits=breaks[2:(length(breaks)-1)],
                       main_title=unit[m], plot_loc=c(0.2,0.8,0.63, 0.77), 
                       title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
    

    
    dev.off()
    
  } #variables
  
  # } #GCM/BC
} #metrics

endCluster()

#stopCluster(cl) #need to put this here, otherwise temporary .grd files craeted by parallel code become unavailable



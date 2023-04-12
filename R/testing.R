library(raster)

cmip=mean(brick(paste0("/g/data/w97/amu561/CABLE_AWRA_comparison/Drought_metrics/", 
                       "CMIP5/historical/mrro/Perc_15/Baseline_1970_2005/Scale_3/CNRM-CM5/r1i1p1/",
                       "CNRM-CM5_r1i1p1_drought_metrics_perc_15_historical_1850_2005.nc"), 
                varname="duration"), na.rm=TRUE)


awra= mean(brick(paste0("/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics/", 
                        "3-month/CCAM/CNRM-CERFACS-CNRM-CM5/drought_metrics_", 
                        "CCAM_CNRM-CERFACS-CNRM-CM5_qtot_rcp45_3.nc"),
                 varname="duration")[[1:200]], na.rm=TRUE)



test=resample(cmip, awra, method="bilinear")

test1=disaggregate(cmip, fac=28)

par(mfcol=c(1,3))
plot(cmip)
plot(test)
plot(test1)



library(ncdf4)
library(raster)


#clear R environment
rm(list=ls(all=TRUE))


#Mengyuan's AWAP gridinfo file with LAI, landsea, patchfrac and iveg removed
gridinfo <- "/g/data/w35/amu561/Steven_CABLE_runs/CABLE_inputs/gridinfo_AWRA_matched_veg_inputs.nc"


nc <- nc_open(gridinfo, write=TRUE)



#######################
### AWRA land cover ###
#######################

tree_frac <- raster(paste0("/g/data/w35/amu561/Steven_CABLE_runs/AWRA_fields/",
                          "AWRA_tree_fraction.nc"))

grass_frac <- 1 - tree_frac


#Create iveg layer that has evergreen broadleaf and C3 grass

iveg <- brick(nl=2, nrows=nrow(tree_frac), ncols=ncol(tree_frac),
              xmn=xmin(tree_frac), xmx=xmax(tree_frac),
              ymn=ymin(tree_frac), ymx=ymax(tree_frac))

iveg[[1]] <- 2
iveg[[2]] <- 6


#Create patchfrac layer
patchfrac <- brick(list(tree_frac, grass_frac))



################
### AWRA LAI ###
################

#AWRA fractional vegetation cover
fcover <- brick(paste0("/g/data/w35/amu561/Steven_CABLE_runs/AWRA_fields/",
                       "AWRA_fractional_vegetation_cover_monthly_climatology_1960_2005.nc"))

#Convert to LAI
#use k=0.7 for shallow-rooted veg and k=0.4 for deep-rooted veg
#following van Dijk 2010 (AWRA technical report)

#CABLE-trunk can't take tiled LAI. Take a weighted mean of k in each pixel
#depending on shallow/deep veg fractions
weighted_k <- calc(patchfrac, function(x) weighted.mean(c(0.4, 0.7), w=x))


lai <- log(1-fcover)/(-weighted_k)

# lai_shallow <- log(1-fcover)/(-0.4)
# 
# lai_deep <- log(1-fcover)/(-0.7)
# 
# 
# #Create LAI layer with both tiles
# lai <- brick(list(lai_deep, lai_shallow))



#####################
### Land-sea mask ###
#####################

landsea <- lai[[1]]
landsea[!is.na(landsea)] <- 0
landsea[is.na(landsea)] <- 1




##########################
### Write new gridinfo ###
##########################


#Define patch dimension (doesn't exist in file)

patchdim <- ncdim_def("patch", units="-", vals=c(1:2))


#Define new variables

lai_var <- ncvar_def("LAI", units="-", dim=list(nc$dim$longitude,
                                                nc$dim$latitude, nc$dim$time), prec="double")

patchfrac_var <- ncvar_def("patchfrac", units="-", dim=list(nc$dim$longitude, nc$dim$latitude, 
                                                            patchdim), prec="double")

iveg_var <- ncvar_def("iveg", units="-", dim=list(nc$dim$longitude, nc$dim$latitude, 
                                                  patchdim), prec="double")

landsea_var <- ncvar_def("landsea", units="-", dim=list(nc$dim$longitude,
                                                          nc$dim$latitude), prec="double")


#Add new variables

nc <- ncvar_add(nc, lai_var)
nc <- ncvar_add(nc, patchfrac_var)
nc <- ncvar_add(nc, landsea_var)
nc <- ncvar_add(nc, iveg_var)


#Write values
ncvar_put(nc, varid=lai_var, vals=as.array(lai))#, start=c(1,1,1), count=c(1,1,12))
ncvar_put(nc, varid=patchfrac_var, vals=as.array(patchfrac))#, start=c(1,1,1), count=c(1,1,2))
ncvar_put(nc, varid=iveg_var, vals=as.matrix(iveg))#, start=c(1,1), count=c(1,1))
ncvar_put(nc, varid=landsea_var, vals=as.matrix(landsea))#, start=c(1,1), count=c(1,1))



nc_close(nc)














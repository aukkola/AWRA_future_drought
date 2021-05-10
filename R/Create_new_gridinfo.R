library(ncdf4)
library(raster)


#clear R environment
rm(list=ls(all=TRUE))

path <- "/g/data/w35/amu561/Steven_CABLE_runs/CABLE_inputs/"


#Mengyuan's AWAP gridinfo file with LAI, landsea, patchfrac and iveg removed
gridinfo <- paste0(path, "/gridinfo_AWRA_matched_veg_inputs.nc")


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

#LAI must be greater than 0.01 or CABLE will go bonkers
#(this threshold is hard-coded into cable_albedo, elsewhere
#is it set by LAI_THRESH which is set to 0.001)
lai[lai <= 0.01] <- 0.011


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

#Mengyuan's file
#landsea <- raster("g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_mask.nc")

landsea <- lai[[1]]
landsea[!is.na(landsea)] <- 0
landsea[is.na(landsea)] <- 1

#Read isoil, it has lots of missing cells. Need to mask these as ocean
isoil <- raster(gridinfo, varname="isoil")

landsea[is.na(isoil)] <- 1



##########################
### Write new gridinfo ###
##########################


#Define patch dimension (doesn't exist in file)

patchdim <- ncdim_def("patch", units="-", vals=c(1:2))

#Adding this as CABLE needs it but it gets removed by CDO when processing
#Mengyuan's file. Normally used for Albedo but not in Mengyuan's file
#(the albedo comes from Mark)
raddim <- ncdim_def("rad", units="-", vals=c(1:3))



#Define new variables

lai_var <- ncvar_def("LAI", units="-", dim=list(nc$dim$longitude,
                                                nc$dim$latitude, nc$dim$time), prec="double")

patchfrac_var <- ncvar_def("patchfrac", units="-", dim=list(nc$dim$longitude, nc$dim$latitude, 
                                                            patchdim), prec="double")

iveg_var <- ncvar_def("iveg", units="-", dim=list(nc$dim$longitude, nc$dim$latitude, 
                                                  patchdim), prec="double")

landsea_var <- ncvar_def("landsea", units="-", dim=list(nc$dim$longitude,
                                                          nc$dim$latitude), prec="double")

#Dummy variable to allow rad dimension to be added to file
dummy_var <- ncvar_def("dummy", units="-", dim=raddim, prec="double")



#Add new variables

nc <- ncvar_add(nc, lai_var)
nc <- ncvar_add(nc, patchfrac_var)
nc <- ncvar_add(nc, landsea_var)
nc <- ncvar_add(nc, iveg_var)
nc <- ncvar_add(nc, dummy_var)


#Write values

#Need to flip and transpose the matrices, otherwise rubbish written into the file. Need to used aperm
#for array and t for matrix
ncvar_put(nc, varid=lai_var, vals=aperm(as.array(flip(lai, direction='y')), c(2,1,3))) #as.array(lai)), 
ncvar_put(nc, varid=patchfrac_var, vals=aperm(as.array(flip(patchfrac, direction='y')), c(2,1,3))) 
ncvar_put(nc, varid=iveg_var, vals=aperm(as.array(flip(iveg, direction='y')), c(2,1,3))) 
ncvar_put(nc, varid=landsea_var, vals=t(as.matrix(flip(landsea, direction='y'))))
ncvar_put(nc, varid=dummy_var, vals=1:3)


nc_close(nc)



# 
# 
# par(mfcol=c(2,2))
# 
# mask <- raster(gridinfo, varname="landsea")
# plot(mask, main="mask")
# 
# lai <- raster(gridinfo, varname="LAI")
# plot(lai, main="LAI")
# 
# iveg <- raster(gridinfo, varname="iveg")
# plot(iveg, main="iveg")
# 
# patch <- raster(gridinfo, varname="patchfrac")
# plot(patch, main="patchfrac")
# 
# 






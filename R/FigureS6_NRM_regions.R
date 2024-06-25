library(raster)


library(RColorBrewer)
library(maptools)
library(maps)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/Steven_CABLE_runs/" 


#Get world shapefile (for masking) and crop to Australia
data("wrld_simpl", package="maptools")
Aus <- crop(wrld_simpl, extent(112.5, 153.75, -43.75, -10))


#NRM regions
nrm_regions <- raster("/g/data/wj02/MISC/NRM_CLUSTERS/NRM_clusters.nc")

#MAsk zero values (ocean)
nrm_regions[nrm_regions ==0] <- NA

#Labels
nrm_labels <- c("Central Slopes", "East Coast", 
                "Murray Basin", "Monsoonal North",
                "Rangelands", "Southern Slopes",
                "S/SW Flatlands",
                "Wet Tropics")


cols <- c("#8dd3c7","#ffffb3","#bebada","#fb8072",
          "#80b1d3","#fdb462","#b3de69","#fccde5")

outdir <- paste0(path, "/Figures/")


### Set up figure ###
png(paste0(outdir, "/FigureS6_NRM_regions.png"),
    height=5, width=5, units="in", res=400)

par(bty="n")

#For some reason, no data for value == 3
breaks = c(0.5, 1.5, 2.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)

#Plot
plot(mask(nrm_regions, Aus), yaxt="n", xaxt="n", col=cols, breaks=breaks, 
     legend=FALSE)

#Australia outline
map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"


#Split legend
cex=0.7
legend(110, -37, legend=nrm_labels[1], fill=cols[1], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(125, -37, legend=nrm_labels[2], fill=cols[2], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(110, -39, legend=nrm_labels[3], fill=cols[3], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(125, -39, legend=nrm_labels[4], fill=cols[4], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(110, -41, legend=nrm_labels[5], fill=cols[5], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(125, -41, legend=nrm_labels[6], fill=cols[6], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(110, -43, legend=nrm_labels[7], fill=cols[7], xpd=NA, horiz=TRUE, cex=cex, bty="n")
legend(125, -43, legend=nrm_labels[8], fill=cols[8], xpd=NA, horiz=TRUE, cex=cex, bty="n")



dev.off()








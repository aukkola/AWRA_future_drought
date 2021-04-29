
new_file="/g/data/w35/amu561/Steven_CABLE_runs/CABLE_inputs/gridinfo_AWRA_matched_veg_inputs.nc"

gridinfo="/g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline/gridinfo_AWAP_OpenLandMap_ELEV_DLCM.nc"

cdo delname,LAI,landsea,patchfrac,iveg $gridinfo $new_file 




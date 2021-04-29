#!/bin/bash

#PBS -m ae
#PBS -P oq98
#PBS -q normal
#PBS -l walltime=01:40:00
#PBS -l mem=3GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35+gdata/er4
#PBS -M a.ukkola@unsw.edu.au


module load cdo


outdir="/g/data/w35/amu561/Steven_CABLE_runs/AWRA_fields"


files=`ls /g/data/er4/awral_nc_versions/from_1950/v61/sim_results/fveg_*.nc`

cdo mergetime $files $outdir"/fveg_allyears.nc"


cdo selyear,1960/2005 $outdir"/fveg_allyears.nc" $outdir"/fveg_1960_2005.nc"

cdo ymonmean -monmean $outdir"/fveg_1960_2005.nc" \
$outdir"/AWRA_fractional_vegetation_cover_monthly_climatology_1960_2005.nc"

rm outdir"/fveg_allyears.nc" $outdir"/fveg_1960_2005.nc"


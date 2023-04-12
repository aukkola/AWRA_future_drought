#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normal
#PBS -l walltime=00:30:00
#PBS -l mem=1GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/oq98+gdata/hh5
#PBS -M a.ukkola@unsw.edu.au

module use /g/data3/hh5/public/modules
module load conda/analysis3-unstable
module load nco


#File path
path="/g/data/w97/amu561/Steven_CABLE_runs/CABLE_outputs"

#Experiments
declare -a experiments=('historical' 'rcp45' 'rcp85')

#CO2/noCO2 directory
declare -a co2_dir=('CO2' 'noCO2')
#co2_dir=`ls $path`

for c in ${co2_dir[@]}
do
  
  gcm_dir=`ls $path/$c`
  
  #Loop through GCMs
  for g in $gcm_dir
  do
    
    #Loop through experiments
    for e in ${experiments[@]}
    do
      
      #Input directory
      indir="$path/${c}/${g}/${e}/r240x120-MRNBC-AWAP/outputs/"
      
      #Find files without path
      files=`ls ${indir}/*.nc` #`find ${indir}/outputs/ -name "*.nc" | sed 's!.*/!!'`   
      
      #Loop through files
      for f in $files
      do
        
        #Input file
        #infile="${indir}/outputs/${f}"
        
        #Fix calendar (change all files, easier than working out leap years...)
        ncatted -h -O -a calendar,time,o,c,"noleap" $f #$infile
        
      done #files
    done #experiments
  done #gcms
done #CO2/noCO2








#CABLE 


ncatted -h -O -a calendar,time,o,c,"noleap" cable_out_2008.nc

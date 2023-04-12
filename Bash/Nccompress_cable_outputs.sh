#!/bin/bash

#PBS -m ae
#PBS -P oq98
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=3GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/oq98+gdata/hh5
#PBS -M a.ukkola@unsw.edu.au

module use /g/data3/hh5/public/modules
module load conda/analysis3-unstable

#File path
path="/g/data/w97/amu561/Steven_CABLE_runs/CABLE_outputs"

#Experiments
declare -a experiments=('historical' 'rcp45' 'rcp85')

declare -a co2_dir=('CO2' 'noCO2')


#CO2/noCO2 directory
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
      indir="$path/${c}/${g}/${e}/r240x120-MRNBC-AWAP/"
      
      #Find files without path
      files=`find ${indir}/outputs/ -name "*.nc" | sed 's!.*/!!'`   
      
      
      #Output directory
      outdir="${path}/${c}_compressed/${g}/${e}/r240x120-MRNBC-AWAP/"
      mkdir -p $outdir $outdir/outputs
      
      
      #Copy restarts, logs and namelists to output location
      cp -r -n $indir/restarts $outdir
      cp -r -n $indir/logs $outdir
      cp -r -n $indir/namelists $outdir
      
      
      #Loop through files
      for f in $files
      do
        
        #Input file
        infile="${indir}/outputs/${f}"
        
        #output file
        outfile="${outdir}/outputs/${f}"      
    
        #Check if exists
        if [ ! -f ${outfile} ]
        then    
        
          #First need to delete froot from input file as CDO can't handle this
          #(CDO needed to check that files have been compressed correctly)
          #NOPE DON'T USE THIS, corrupts the file. Very annoying and not sure why
          #cdo delvar,froot $infile $infile
        
          #Compress file (can't use the -p flag as the cdo command that it calls fails
          #because CDO doesn't like the froot variable in the files)
          nccompress $infile
          
          #Move file to output location
          mv ${indir}/outputs/tmp.nc_compress/${f} $outfile
          
        fi
        
      done #files
    done #experiments
  done #gcms
done #CO2/noCO2










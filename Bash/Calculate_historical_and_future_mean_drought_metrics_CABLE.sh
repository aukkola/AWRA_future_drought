#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normal
#PBS -l walltime=02:40:00
#PBS -l mem=3GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97
#PBS -M a.ukkola@unsw.edu.au




#Set path
path="/g/data/w97/amu561/Steven_CABLE_runs/"

#Find scales
#scale=`ls $path/drought_metrics`
declare -a scale=("3" "24")


#Set experiments
declare -a exp=("historical" "rcp45")

#Set variables
declare -a vars=("qtot" "sm")

#Set metrics
declare -a metrics=("duration" "rel_intensity")


#Loop through scales
for s in ${scale[@]}
do
  
  #Find bias correction methods
  co2_options=`ls "${path}/drought_metrics_CABLE/${s}-month"`
  
  #Loop through CO2 optionss
  for b in $co2_options
  do
    
    #Find GCMs
    gcms=`ls "${path}/drought_metrics_CABLE/${s}-month/${b}"`

    #Loop through GCMs
    for g in $gcms
    do

      #Loop through variables
      for v in ${vars[@]}
      do
        
        #Loop through experiments
        for e in ${exp[@]}
        do
          
          #Set time period
          if [ $e = "historical" ]; then
            start_yr=1970
            end_yr=2005
          else
            start_yr=2064
            end_yr=2099
          fi
          
          
          echo "exp: ${e}, start yr: ${start_yr}, end yr: ${end_yr}"
          
          #Find data file
          inpath="${path}/drought_metrics_CABLE/${s}-month/${b}/${g}/"
          infile="${inpath}/drought_metrics_CABLE_${b}_${g}_${v}.nc"
  
          
          #Loop through metrics
          for m in ${metrics[@]}
          do
            
            outdir="${path}/Mean_drought_metrics_CABLE/scale_${s}/${v}/${b}/${g}/"
            mkdir -p $outdir
            
            outfile="${outdir}/Mean_CABLE_${m}_${v}_${b}_${g}_${e}_${start_yr}_${end_yr}_scale_${s}.nc"
            
            if [ ! -f ${outfile} ]; then
              #Select metric and time period and calculate mean
              cdo -L timmean -selyear,$start_yr/$end_yr -selvar,$m $infile $outfile 
          
            fi
          
          done #metrics
        
        done #experiments
    
      done #vars
    
    done #GCMS
  done #bc_methods

done #scale

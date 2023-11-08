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
declare -a exp=("historical" "rcp45" "rcp85")

#Set variables
declare -a vars=("pr" "qtot" "sm")

#Set metrics
declare -a metrics=("duration" "rel_intensity" "frequency")


#Loop through scales
for s in ${scale[@]}
do
  
  #Find bias correction methods
  bc_methods=`ls "${path}/drought_metrics/${s}-month"`
  
  #Loop through bc methods
  for b in $bc_methods
  do
    
    #Find GCMs
    gcms=`ls "${path}/drought_metrics/${s}-month/${b}"`

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
            pattern="rcp45" #file pattern to read, use rcp for historical as don't have a separate file for it
          else
            start_yr=2064
            end_yr=2099
            pattern=$e
          fi
          
          
          echo "exp: ${e}, start yr: ${start_yr}, end yr: ${end_yr}"
          
          #Find data file
          inpath="${path}/drought_metrics/${s}-month/${b}/${g}/"
          infile="${inpath}/drought_metrics_${b}_${g}_${v}_${pattern}_${s}.nc"
  
          
          #Loop through metrics
          for m in ${metrics[@]}
          do
            
            outdir="${path}/Mean_drought_metrics/scale_${s}/${e}/${v}/${b}/${g}/"
            mkdir -p $outdir
            
            outfile="${outdir}/Mean_${m}_${v}_${b}_${g}_${e}_${start_yr}_${end_yr}_scale_${s}.nc"
            
            if [ ! -f ${outfile} ]; then
              
              if [ $m = "frequency" ]; then
                
                no_yrs="$(($end_yr-$start_yr+1))"
                
                #Select metric and time period and calculate fraction of time under drought
                #sum timing (i.e. binary 0/1 drought/no-drought), then multiply the number
                #of years and months and finally convert from fraction to percentage by multiplying by 100
                cdo -L expr,"timing=timing/($no_yrs*12)*100" -timsum -selyear,$start_yr/$end_yr -selvar,timing $infile $outfile 
          
              else
                #Select metric and time period and calculate mean
                cdo -L timmean -selyear,$start_yr/$end_yr -selvar,$m $infile $outfile 
            
              fi
          
            fi
          
          done #metrics
        
        done #experiments
    
      done #vars
    
    done #GCMS
  done #bc_methods

done #scale

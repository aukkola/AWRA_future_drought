#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normal
#PBS -l walltime=04:40:00
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
declare -a metrics=("frequency" "rel_intensity_by_month")


#Loop through scales
for s in ${scale[@]}
do
  
  echo "scale ${s}"
  
  #Find bias correction methods
  co2_options=`ls "${path}/drought_metrics_CABLE/${s}-month"`
  
  #Loop through bc methods
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
          
          echo "exp: ${e}, start yr: ${start_yr}, end yr: ${end_yr}"

          #Set time period
          if [ $e = "historical" ]; then
            start_yr=1970
            end_yr=2005
          else
            start_yr=2064
            end_yr=2099
          fi
                    
          
          #Find data file
          inpath="${path}/drought_metrics_CABLE/${s}-month/${b}/${g}/"
          infile="${inpath}/drought_metrics_CABLE_${b}_${g}_${v}.nc"
  
          
          #Loop through metrics
          for m in ${metrics[@]}
          do
            
            echo "metric: ${m}, variable: ${v}, CO2: ${b}, gcm: ${g}"


            outdir="${path}/Seasonal_drought_metrics_CABLE/scale_${s}/${e}/${v}/${b}/${g}/"
            mkdir -p $outdir
            
            outfile="${outdir}/Mean_${m}_${v}_${b}_${g}_${start_yr}_${end_yr}_scale_${s}"


            #Need to take time-sum for frequency to work out time under drought
            #Take mean for other variables

            if [ $m = "frequency" ]; then
              
              #Calculate separately for each season
              #use timing variable

              #Check if exists (only do summer file for simplicity...)
              if [ ! -f ${outfile}"_summer.nc" ]
              then

                #Summer
                cdo -L timsum -selseas,DJF -selyear,$start_yr/$end_yr -selvar,timing $infile ${outfile}"_summer.nc"
                #Autumn
                cdo -L timsum -selseas,MAM -selyear,$start_yr/$end_yr -selvar,timing $infile ${outfile}"_autumn.nc"
                #Winter
                cdo -L timsum -selseas,JJA -selyear,$start_yr/$end_yr -selvar,timing $infile ${outfile}"_winter.nc"
                #Spring
                cdo -L timsum -selseas,SON -selyear,$start_yr/$end_yr -selvar,timing $infile ${outfile}"_spring.nc"
                
              fi
              
            else
              
              if [ ! -f ${outfile}"_summer.nc" ]
              then

                #Calculate separately for each season

                #Summer
                cdo -L timmean -selseas,DJF -selyear,$start_yr/$end_yr -selvar,$m $infile ${outfile}"_summer.nc"
                #Autumn
                cdo -L timmean -selseas,MAM -selyear,$start_yr/$end_yr -selvar,$m $infile ${outfile}"_autumn.nc"
                #Winter
                cdo -L timmean -selseas,JJA -selyear,$start_yr/$end_yr -selvar,$m $infile ${outfile}"_winter.nc"
                #Spring
                cdo -L timmean -selseas,SON -selyear,$start_yr/$end_yr -selvar,$m $infile ${outfile}"_spring.nc"
              fi
            fi

                          
          done #metrics
        
        done #experiments
    
      done #vars
    
    done #GCMS
  done #bc_methods

done #scale

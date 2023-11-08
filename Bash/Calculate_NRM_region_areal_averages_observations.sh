#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normal
#PBS -l walltime=5:40:00
#PBS -l mem=3GB
#PBS -l ncpus=1
#PBS -l jobfs=6Gb
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97
#PBS -M a.ukkola@unsw.edu.au

module load cdo
module load R

#Set path
path="/g/data/w97/amu561/Steven_CABLE_runs"

#Set scale
scale=3

#Set metrics
declare -a metrics=("duration" "rel_intensity" "timing")

#Set variables
declare -a vars=("pr" "qtot" "sm")


#9 NRM regions in total (numbered 1-9)

#NRM regions
nrm_file="$path/NRM_means/NRM_regions_masked.nc"


#Loop through variables
for v in ${vars[@]}
do

  #Select AGCD data for pr and AWRA reference run for qtot and sm
  if [ $v == "pr" ]
  then
    infile="${path}/drought_metrics_AGCD/${scale}-month/drought_metrics_AGCD_precip_1900_2021_baseline_1970_2005_scale_${scale}.nc"
  else
    infile="${path}/drought_metrics_AWRA_ref/${scale}-month/drought_metrics_AWRA_ref_${v}_scale_${scale}_1960_2020.nc"
  fi
  

  #Loop through metrics
  for m in ${metrics[@]}
  do
    
    #output directory
    outdir=$path"/NRM_means/data_obs/scale_"$scale"/"$m"/"$v"/"
    mkdir -p $outdir
    
    #Temp file
    temp_file="$outdir/temp.nc"

    #Select metric

    if [ $v == "pr" ]
    then
    
    cdo -L invertlat -selvar,$m $infile $temp_file

    echo "temp_file $temp_file"
    
    
    cat > Rscript.R << EOF 
    library(raster)

    agcd <- brick("${temp_file}")
    nrm  <- raster("${nrm_file}")

    data <- crop(agcd, nrm)

    writeRaster(data, "${temp_file}", varname="$m", overwrite=TRUE)

EOF

    Rscript Rscript.R
    
    else
    cdo -L -selvar,$m $infile $temp_file

    fi


    #Loop through NRM regions
    for nrm in $(seq 1 9)
    do
      
      temp_masked="$outdir/temp_masked.nc"
      
      #Mask NRM region
      cdo -L div $temp_file -eqc,$nrm $nrm_file $temp_masked

      #Output file
      outfile="${outdir}/Areal_mean_${m}_${v}_scale_${scale}_region${nrm}.nc"

      #Check if file exists, skip otherwise
      if [ ! -f $outfile ]
      then
    
        #If variable is timing, take areal sum
        if [ $m == "timing" ]
        then
          cdo -L fldsum $temp_masked $outfile
        #Else take the mean
        else
          cdo -L fldmean $temp_masked $outfile
        fi
        
      fi    
      
      
      #Remove temp file
      rm $temp_masked
  
  
    done #NRM regions

    rm $temp_file 
  
  done #metric
done #variable

  

  
  

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

module load cdo

#Set path
path="/g/data/w97/amu561/Steven_CABLE_runs"

#Set scale
scale=3

#Set metrics
declare -a metrics=("duration" "rel_intensity" "timing")

#Set variables
declare -a vars=("pr" "qtot" "sm")

#Set variables
declare -a scenarios=("rcp45" "rcp85")


#9 NRM regions in total (numbered 1-9)

#NRM regions
nrm_file="$path/NRM_means/NRM_regions_masked.nc"


#Find BC methods
bc_methods=`ls "${path}/drought_metrics/${scale}-month"`



#Loop through scenarios
for scen in ${scenarios[@]}
do

  #Loop through bc methods
  for bc in $bc_methods
  do
    
    #Find GCMs
    gcms=`ls "${path}/drought_metrics/${scale}-month/${bc}"`
    
    #Loop through GCMs
    for gcm in $gcms
    do
        #Loop through variables
        for v in ${vars[@]}
        do

          infile=`ls ${path}/drought_metrics/${scale}-month/${bc}/${gcm}/*_${v}_${scen}_${scale}.nc`
    
          #Loop through metrics
          for m in ${metrics[@]}
          do
            
            #output directory
            outdir=$path"/NRM_means/data/scale_"$scale"/"$m"/"$scen"/"$v"/"$bc"/"$gcm
            mkdir -p $outdir
            
            #Temp file
            temp_file="$outdir/temp.nc"

            #Select metric
            cdo -L invertlat -selvar,$m $infile $temp_file
      
      
            #Loop through NRM regions
            for nrm in $(seq 1 9)
            do
              
              temp_masked="$outdir/temp_masked.nc"
              
              #Mask NRM region
              cdo -L div $temp_file -eqc,$nrm $nrm_file $temp_masked
       
              #Output file
              outfile="${outdir}/Areal_mean_${m}_${bc}_${gcm}_${scen}_${v}_scale_${scale}_region${nrm}.nc"
       
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
    done #gcm
  done #bc method
done #scenario

  

  
  

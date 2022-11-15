#!/bin/bash

#PBS -m ae
#PBS -P oq98
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=60GB
#PBS -l jobfs=3Gb
#PBS -l ncpus=48
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/wd9+scratch/w97+gdata/wj02
#PBS -M a.ukkola@unsw.edu.au

module unload openmpi
module add intel-mpi/2019.5.281
module add netcdf/4.7.4p


### Step 3 wall time ###

#Maximum wall time is 48hrs, do not change it

#(3.5 hours per year)


#################
### Run CABLE ###
#################


#Met input directory
met_indir=$wg_out_path"/"${model}"/"${experiment}"/"${bc_method}"/"


#Change to CABLE directory and compile
cd $cable_src_path'/trunk_31Mar2021/offline'

startyr=$year


#Do in batches of 5 years
while [ $year -le $((startyr+4)) -a $year -le $EndYear ]
do
  

  echo "Step 3: Running CABLE (increasing CO2) for $year #-----------------------"

  #Set previous year
  prev_year=$((year-1))
  
  
  #Set output files
  logfile=$logs_CO2"/cable_log_${year}.txt"
  outfile=$outs_CO2"/cable_out_${year}.nc"
  restart_in=$restarts_CO2"/restart_${prev_year}.nc"
  restart_out=$restarts_CO2"/restart_${year}.nc"
  namelist=$nml_CO2"/cable_${year}.nml"


  #Get CO2 concentration for the year (second command removes a line ending character)
  co2=`echo "${co2_file[$year]}" | tr '\r' ' ' `

  #Create namelist
  sh ./create_cable-nml_co2.sh -y $year -l $logfile -o $outfile -i $restart_in -r $restart_out -c $co2 -m $met_indir

  #Run CABLE
  mpirun -n 48 ./cable-mpi ./cable_on.nml

  #Copy namelist
  cp cable_on.nml $namelist
  
  year=$((year+1))
  
done


#NEED TO TEST FOR CO2 and noCO2 before removing !!!!!!!!!!!!!!!!!!!!!

#!!!!!

#If not end year, qsub again
if [ $year_ind -lt $EndYr ]; then
  qsub 


#If end year, tidy up
else 

  ### Tidy up ###
  
  echo "Step 5: Tidying up #-----------------------"


  #Check that have same number of CABLE output files as
  #forcing files
  no_years=`seq $startYr $endYr | wc -l`
  cable_outs_CO2=`ls $outs_CO2"/*.nc" | wc -l`

  #If they don't match, stop without deleting forcing files

  if [ $no_years -ne $cable_outs_CO2 ]; then 
    echo "ERROR: The number of CABLE output files (inc. CO2) does not match the number of years"
    exit 1
  fi


  cable_outs_noCO2=`ls $outs_noCO2"/*.nc" | wc -l`

  if [ $no_years -ne $cable_outs_noCO2 ]; then 
    echo "ERROR: The number of CABLE output files (const. CO2) does not match the number of years"
    exit 1
  fi

  #Remove WG forcing files
  rm -r $met_indir

fi










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
#PBS -l storage=gdata/w97+gdata/wd9+scratch/w97+gdata/wj02+gdata/oq98
#PBS -M a.ukkola@unsw.edu.au

module unload openmpi
module unload netcdf
module add intel-mpi/2019.5.281
module add netcdf/4.6.3


### Step 3 wall time ###

#Maximum wall time is 48hrs, do not change it

#(3.5 hours per year)



################
### Settings ###
################

#Grab path where run script is
rundir=`pwd`

cd $rundir 

#CABLE output path (where CABLE outputs will be stored)
cable_out_path_CO2=$path"/CABLE_outputs/CO2/"$model/$experiment/$bc_method
cable_out_path_noCO2=$path"/CABLE_outputs/noCO2/"$model/$experiment/$bc_method



#Set CO2 file
if [[ $experiment == "historical" || $experiment == "rcp45" ]]; then 
  co2_file=(`cat "/g/data/w97/amu561/Steven_CABLE_runs/CO2_concentrations/CABLE_Annual_CO2_concentration_1_2100_RCP4.5.csv"`)

elif [[ $experiment == "rcp85" ]]; then
  co2_file=(`cat "/g/data/w97/amu561/Steven_CABLE_runs/CO2_concentrations/CABLE_Annual_CO2_concentration_1_2100_RCP8.5.csv"`)

fi


#################
### Run CABLE ###
#################


#Create output folders (separately for increasing and constant CO2)
mkdir -p $cable_out_path_CO2
mkdir -p $cable_out_path_noCO2


### Create subdirectories ###

#CO2
outs_CO2=$cable_out_path_CO2"/outputs/"
restarts_CO2=$cable_out_path_CO2"/restarts/"
nml_CO2=$cable_out_path_CO2"/namelists/"
logs_CO2=$cable_out_path_CO2"/logs/"

mkdir $outs_CO2
mkdir $restarts_CO2
mkdir $nml_CO2
mkdir $logs_CO2

#No CO2
outs_noCO2=$cable_out_path_noCO2"/outputs/"
restarts_noCO2=$cable_out_path_noCO2"/restarts/"
nml_noCO2=$cable_out_path_noCO2"/namelists/"
logs_noCO2=$cable_out_path_noCO2"/logs/"

mkdir $outs_noCO2
mkdir $restarts_noCO2
mkdir $nml_noCO2
mkdir $logs_noCO2


#Met input directory
met_indir=$wg_out_path"/"${model}"/"${experiment}"/"${bc_method}"/"

#First year of this iteration, used to only run 5yrs at a time below
first_yr=$year

#Run first 5 years
while [ $year -le $((first_yr+10)) -a $year -le $endYr ]
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

  #Check that restart file exists (except for start year of historical expt)
  #stop if not to avoid cable running without a restart
  
  #If restart doesn't exist
  if [ ! -f $restart_in ]
  then
    if [ $experiment == "historical" -a $year -eq $((startYr)) ]
    then
      echo "restart doesn't exist, first year of historical experiment"
    else
      echo "ERROR: restart_in does not exist"
      exit 1
    fi
  fi

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


#Then submit next 5 years
#line break problems again, having it all on one line...
cd $rundir

qsub -v "path=$path","wg_out_path=$wg_out_path","model=$model","experiment=$experiment",\
"bc_method=$bc_method","startYr=$startYr","endYr=$endYr","year=$year","cable_src_path=$cable_src_path" Run_cable_step3.sh 




# #Loop through years
# for year in $(seq $startYr $endYr)
# do
# 
#   echo "Step 3: Running CABLE (increasing CO2) for $year #-----------------------"
# 
#   #Met input directory
#   met_indir=$wg_out_path"/"${model}"/"${experiment}"/"${bc_method}"/"
# 
# 
# 
#   ######################
#   ### Increasing CO2 ###
#   ######################
# 
#   #Set output files
#   logfile=$logs_CO2"/cable_log_${year}.txt"
#   outfile=$outs_CO2"/cable_out_${year}.nc"
#   restart_in=$restarts_CO2"/restart_${prev_year}.nc"
#   restart_out=$restarts_CO2"/restart_${year}.nc"
#   namelist=$nml_CO2"/cable_${year}_${gw_tag}.nml"
# 
# 
#   #Get CO2 concentration for the year (second command removes a line ending character)
#   co2=`echo "${co2_file[$year]}" | tr '\r' ' ' `
# 
#   #Create namelist
#   sh ./create_cable-nml_co2.sh -y $year -l $logfile -o $outfile -i $restart_in -r $restart_out -c $co2 -m $met_indir
# 
#   #Run CABLE
#   mpirun -n 48 ./cable-mpi ./cable_on.nml
# 
#   #Copy namelist
#   cp cable_on.nml $namelist
# 
# done
# 


# #Loop through years
# for year in $(seq $startYr $endYr)
# do
# 
#   echo "Step 4: Running CABLE (constant CO2) for $year #-----------------------"
# 
# 
#   ####################
#   ### Constant CO2 ###
#   ####################
# 
#   #Run with CO2 set to 1960 level
# 
#   #Set output files
#   logfile=$logs_noCO2"/cable_log_${year}.txt"
#   outfile=$outs_noCO2"/cable_out_${year}.nc"
#   restart_in=$restarts_noCO2"/restart_${prev_year}.nc"
#   restart_out=$restarts_noCO2"/restart_${year}.nc"
#   namelist=$nml_noCO2"/cable_${year}.nml"
# 
# 
#   #Get CO2 concentration for 1960
#   co2=`echo "${co2_file[1960]}" | tr '\r' ' ' `
# 
#   #Create namelist
#   sh ./create_cable-nml_co2.sh -y $year -l $logfile -o $outfile -i $restart_in -r $restart_out -c $co2 -m $met_indir
# 
#   #Run CABLE
#   mpirun -n 48 ./cable-mpi ./cable_on.nml
# 
#   prev_year=$year
# 
#   #Copy namelist
#   cp cable_on.nml $namelist
# 
# done


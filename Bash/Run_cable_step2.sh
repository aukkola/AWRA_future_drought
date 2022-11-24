#!/bin/bash

#PBS -m ae
#PBS -P oq98
#PBS -q normal
#PBS -l walltime=15:00:00
#PBS -l mem=30GB
#PBS -l jobfs=1Gb
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+scratch/w97
#PBS -M a.ukkola@unsw.edu.au


module unload netcdf
module unload openmpi
module add intel-mpi/2019.5.281
module add netcdf/4.6.3


### Step 2 wall time ###

#historical run: 
#walltime=15:00:00

#future run:
#walltime=40:00:00


#(12 min per year)



################
### Settings ###
################

#Set start and end years
if [[ $experiment == "historical" ]]; then 
  startYr=1960
  endYr=2005

elif [[ $experiment == "rcp45" || $experiment == "rcp85" ]]; then
  startYr=2006
  endYr=2099

fi

#Grab path where run script is
rundir=`pwd`


#############################
### Run weather generator ###
#############################

echo "Step 2: Running weather generator #-----------------------"

cd $wg_path

#Pass model, experiment, bc_method, in and out paths through command line
#(in this order)

#Input path (don't change, needs to match what's set in R code)
wg_in_path=$scratch_path"/CABLE_inputs/Weather_generator_inputs/"

#Output path
wg_out_path=$scratch_path"/CABLE_inputs/Weather_generator_outputs/"


#Run weather generator
echo "Running WG"
./awap_to_netcdf $model $experiment $bc_method $wg_in_path $wg_out_path


#Check that produces the correct number of files
#were produced by weather generator
check_path=$wg_out_path"/"${model}"/"${experiment}"/"${bc_method}"/Rainf"

no_years=`seq $startYr $endYr | wc -l`
wg_outs=`ls $check_path/*.nc | wc -l`

#If they don't match, stop without deleting forcing files

if [ $no_years -ne $wg_outs ]; then 
  echo "ERROR: The number of WG output files does not match the number of years"
  exit 1
fi


#Remove daily .flt files
rm -r $wg_in_path"/"${model}"/"${experiment}"/"${bc_method}


echo "Step 2: Finished running weather generator #-----------------------"



#################################
### Submit CABLE run (step 3) ###
#################################

echo "Submitting CABLE run job #-----------------------"

cd $rundir


qsub -v "path=$path","wg_out_path=$wg_out_path","model=$model","experiment=$experiment",\
"bc_method=$bc_method","startYr=$startYr","endYr=$endYr","year=$startYr",\
"cable_src_path=$cable_src_path" Run_cable_step3.sh 







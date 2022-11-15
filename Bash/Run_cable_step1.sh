#!/bin/bash

#PBS -m ae
#PBS -P oq98
#PBS -q normal
#PBS -l walltime=1:50:00
#PBS -l mem=15GB
#PBS -l jobfs=100Mb
#PBS -l ncpus=6
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/wd9+scratch/w97+gdata/wj02
#PBS -M a.ukkola@unsw.edu.au


module unload netcdf
module unload openmpi
module add intel-mpi/2019.5.281
module add netcdf/4.6.3
module load gdal
module load R


### THINGS TO DO ###

# - set wall time for each step
# - change model, experiment and bc_method in this file

#------------------------------------------------------


### Step 1 wall time ###

#historical run: 
#walltime=1:50:00
#mem=15GB

#future run:
#walltime=3:30:00
#mem=25GB

#(2 min per year)


################
### Settings ###
################

#CHANGE model, experiment and bc method options here

model="CNRM-CERFACS-CNRM-CM5" #CNRM-CERFACS-CNRM-CM5,  CSIRO-BOM-ACCESS1-0, MIROC-MIROC5, NOAA-GFDL-GFDL-ESM2M
experiment="historical"   #historical, rcp45, rcp85
bc_method="r240x120-MRNBC-AWAP"

#Set base path
path="/g/data/w97/$USER/Steven_CABLE_runs"




#---- Don't need to change past this if using my file structure

#Scratch path (where temp files will be stored)
scratch_path="/scratch/w97/$USER/Steven_CABLE_runs" 
mkdir -p $scratch_path

#Paths to weather generator and CABLE code
#(load weather generator and CABLE codes to these locations)
wg_path=$path"/AWAP_to_netcdf/"
cable_src_path=$path"/CABLE_source/"


#Path to AWRA inputs (don't change)
awra_path="/g/data/wj02/COMPLIANT/HMINPUT/output/AUS-5/BoM/"



#Grab path where run script is
rundir=`pwd`



###########################################
### Compile weather generator and CABLE ### 
###########################################

cd $wg_path

#Compile code
echo "Compiling"
./wrapper.sh


#Change to CABLE directory and compile
cd $cable_src_path'/trunk_31Mar2021/offline'

./build_mpi.ksh




###########################################
### Create inputs for weather generator ###
###########################################

#All AWRA years in one file so can't run in year loop

echo "Step 1: Creating inputs for weather generator #-----------------------"


#Process AWRA inputs for weather generator
Rscript $path"/scripts/R/Convert_netcdf_to_flt.R" $path $awra_path $scratch_path \
$model $experiment $bc_method


echo "Finished creating inputs for weather generator #-----------------------"



#############################################
### Submit weather generator job (step 2) ###
#############################################

echo "Submitting weather generator job #-----------------------"

#cd $path"/scripts/Bash/"
#cd to run directory to submit next script
cd $rundir

#copy executable and namelist here
cp $cable_src_path/trunk_31Mar2021/offline/cable-mpi .
cp $cable_src_path/trunk_31Mar2021/offline/create_cable-nml_co2.sh .



qsub -v "path=${path}","scratch_path=${scratch_path}","model=${model}","wg_path=${wg_path}",\
"experiment=${experiment}","bc_method=${bc_method}","cable_src_path=${cable_src_path}" Run_cable_step2.sh 



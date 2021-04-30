#!/bin/bash

#PBS -m ae
#PBS -P oq98
#PBS -q normal
#PBS -l walltime=25:00:00
#PBS -l mem=60GB
#PBS -l jobfs=3Gb
#PBS -l ncpus=48
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35+gdata/wd9+scratch/w35+gdata/wj02
#PBS -M a.ukkola@unsw.edu.au


module unload netcdf
module unload openmpi
module add intel-mpi/2019.5.281
module add netcdf/4.6.3
module load gdal
module load R



################
### Settings ###
################

user="amu561"

model="CNRM-CERFACS-CNRM-CM5"
experiment="historical"   #one of historical, rcp45, rcp85
bc_method="CSIRO-CCAM-r3355-r240x120-ISIMIP2b-AWAP"


#Path to scripts (load scripts to this location)
path="/g/data/w35/$user/Steven_CABLE_runs"


#Scratch path (where temp files will be stored)
scratch_path="/scratch/w35/$user/Steven_CABLE_runs" 


#Paths to weather generator and CABLE code
#(load weather generator and CABLE codes to these locations)
wg_path=$path"/AWAP_to_netcdf/"
cable_src_path=$path"/CABLE_source/"


#CABLE output path (where CABLE outputs will be stored)
cable_out_path=$path"/CABLE_outputs/$model/$experiment/$bc_method"


#Path to AWRA inputs (don't change)
awra_path="/g/data/wj02/COMPLIANT/HMINPUT/output/AUS-5/BoM/"



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


#Create output folder
mkdir -p $cable_out_path


#Create subdirectories
outs=$cable_out_path"/outputs/"
restarts=$cable_out_path"/restarts/"
nml=$cable_out_path"/namelists/"
logs=$cable_out_path"/logs/"

mkdir $outs
mkdir $restarts
mkdir $nml
mkdir $logs


#Set CO2 file
if [[ $experiment == "historical" || $experiment == "rcp45" ]]; then 
  co2_file=(`cat "/g/data/w35/amu561/Steven_CABLE_runs/CO2_concentrations/CABLE_Annual_CO2_concentration_1_2100_RCP4.5.csv"`)

elif [[ $experiment == "rcp85" ]]; then
  co2_file=(`cat "/g/data/w35/amu561/Steven_CABLE_runs/CO2_concentrations/CABLE_Annual_CO2_concentration_1_2100_RCP8.5.csv"`)

fi


#Set start and end years
if [[ $experiment == "historical" ]]; then 
  startYr=1960
  endYr=2005

elif [[ $experiment == "rcp45" || $experiment == "rcp85" ]]; then
  startYr=2006
  endYr=2099

fi



###########################################
### Create inputs for weather generator ###
###########################################

#All AWRA years in one file so can't run in year loop

echo "Step 1: Creating inputs for weather generator #-----------------------"


#Process AWRA inputs for weather generator
Rscript $path"/scripts/R/Convert_netcdf_to_flt.R" $path $awra_path $scratch_path \
$model $experiment $bc_method


echo "Finished creating inputs for weather generator #-----------------------"



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




#################
### Run CABLE ###
#################

#Loop through years
for year in $(seq $startYr $endYr)
do
  

  echo "Step 3: Running CABLE for $year #-----------------------"


  #Change to CABLE directory and compile
  cd $cable_src_path'/trunk_31Mar2021/offline'
  

  #Set output files
  logfile=$logs"/cable_log_${year}.txt"
  outfile=$outs"/cable_out_${year}.nc"
  restart_in=$restarts"/restart_${prev_year}.nc"
  restart_out=$restarts"/restart_${year}.nc"
  namelist=$nml"/cable_${year}_${gw_tag}.nml"



  #Need to pass PSurf file depending on whether leap year or not

  #Get CO2 concentration for the year
  co2=`echo "${co2_file[$year]}"`


  #Met input directory
  met_indir=$wg_out_path"/"${model}"/"${experiment}"/"${bc_method}"/"

  #Create namelist
  sh ./create_cable-nml_co2.sh -y $year -l $logfile -o $outfile -i $restart_in -r $restart_out -c $co2 -m $met_indir


  #Run CABLE
  mpirun -n 48 ./cable-mpi ./cable_on.nml

  

  ###############
  ### Tidy up ###
  ###############

  echo "Step 4: Tidying up #-----------------------"


done


#Check that have same number of CABLE output files as
#forcing files
no_years=`seq $startYr $endYr | wc -l`
cable_outs=`ls $outs"/*.nc" | wc -l`

#If they don't match, stop without deleting forcing files

if [ $no_years -ne $cable_outs ]; then 
  echo "ERROR: The number of output files does not match the number of years"
  exit 1
fi

# cable_outs=`ls $outs"/*${year}*.nc" | wc -l`
# 
# if [ $cable_outs -ne 1 ]; then 
#   echo "ERROR: The number of output files does not match the number of years"
#   exit 1
# fi
# 



#Remove WG forcing files
rm -r $wg_out_path

#Remove daily .flt files
rm -r $wg_in_path











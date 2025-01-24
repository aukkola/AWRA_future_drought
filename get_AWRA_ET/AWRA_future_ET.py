## Adapted code from Anna Ukkola 

#################################
### IMPORT NECESSARY PACKAGES ###
#################################

from netCDF4 import Dataset,num2date 
import numpy as np
import glob
import sys 
import os
import datetime
import xarray as xr

#################
### Set paths ###
#################

lib_path  = "/g/data/w97/amu561/Steven_CABLE_runs/scripts/drought_scripts/functions"

##### ALTERTED #####
#data_path = "/g/data/w97/amu561/Steven_CABLE_runs/"
data_path = "/g/data/wj02/COMPLIANT_PUBLISHED/"

scratch_path = '/scratch/w97/amu561/'


# Add lib_path to os directory
sys.path.append(os.path.abspath(lib_path))

#import all droughtmetric calcuation functions
from drought_metrics import *

########################
### DEFINE VARIABLES ###
########################

### Bias Correction ###
bias_corr=str(sys.argv[1])

### Set Models ###
model=str(sys.argv[2])

### Set scenarios ###
scenario=str(sys.argv[3])

### Set variable ###
variable=str(sys.argv[4])

### Set drought metric conditions ###
return_all_tsteps=True

### Set percentile for drought threshold ###
### Set scale for month aggregation ###
scale=3
##########################################
##########################################
##########################################
##########################################


### Set to monthly threshold calc ###
monthly=True

### Set additional reference data to false ###
obs_ref = False
obs_var = variable  

### Set historical refernce period ###
baseline=[1970,2005]

##########################
### FILE PREPROCESSING ###
##########################

### Define file locations ###
output_variable=['etot', 'qtot', 's0', 'ss'] #['Qsb','SoilMoist']
input_variable=['pr']

#compliant=['MRNBC']
#non_compliant=['RAW-GCM','NOBC-CCAM']

# Any output data from CCAM has a slightly differenet path
ccam_add=""
if bias_corr=="CCAM":
    ccam_add="CSIRO-CCAM-r3355-"

awra_add=""
if variable in output_variable:
    awra_add='AWRALv6-1-'

if variable in input_variable: #and bias_corr in compliant:
    data_path_var= data_path + "/HMINPUT/output/AUS-5/BoM"

if variable in output_variable: #and bias_corr in compliant:
    data_path_var= data_path + "/HMOUTPUT/output/AUS-5/BoM"

### Get historical and future simulations ###

#historical
if bias_corr=="CCAM":
    files_1_string=str(data_path_var + '/' + awra_add + model + '/historical/r1i1p1/' +
                       ccam_add + 'r240x120-ISIMIP2b-AWAP/latest/day/' +
                       variable+ '/*1960*.nc')
else:
    files_1_string=str(data_path_var + '/' + awra_add + model +'/historical/r1i1p1/' +
                       ccam_add + 'r240x120-' + bias_corr + '-AWAP/latest/day/' +
                       variable + '/*1960*.nc')

#Files to merge                    
files_to_merge1=glob.glob(files_1_string)

#Future
if bias_corr=="CCAM":
    files_2_string=str(data_path_var + '/' + awra_add + model + '/' + scenario + 
                    '/r1i1p1/' + ccam_add + 'r240x120-ISIMIP2b-AWAP/' + 
                    'latest/day/' + variable + '/*2006*.nc')
else:
    files_2_string=str(data_path_var + '/' + awra_add + model + '/' + scenario + 
                       '/r1i1p1/' + ccam_add + 'r240x120-' + bias_corr + 
                       '-AWAP/latest/day/' + variable + '/*2006*.nc')

files_to_merge2=glob.glob(files_2_string)

files_to_merge1.extend(files_to_merge2)


### Create temporary directory ###

#daily temporary file
temp_dir_path = f"/scratch/w97/amu561/temp/{bias_corr}{model}{scenario}{variable}{str(scale)}"
os.system("mkdir -p " + temp_dir_path)

#monthly temporary file
mon_temp_path = str(scratch_path + "/monthly_sums_AWRA_GCM/")
os.system("mkdir -p " + mon_temp_path)


### Location of output file ###
files= str(mon_temp_path + "/" + bias_corr + "_" + 
           model + "_" + scenario + "_" + variable + ".nc")


#Some duplicate time steps, need to skip these when merging
os.system("export SKIP_SAME_TIME=1")

#The code should be rewritten, it's a bit messy. But in the meantime,
#need to define averaging method separately for soil moisture (as need to take
#monthly mean rather than sum)
if variable == "s0" or variable == "ss":
    fun="monmean"
else: 
    fun="monsum"


#Output file    
outfile = str(scratch_path + "/monthly_sums/" + bias_corr + "_" + model + "_" + 
              scenario + "_" + variable + ".nc")

#If output file doesn't alreade exist, create it. Else skip 
if not os.path.isfile(outfile):

    for file_ms in range(len(files_to_merge1)):
    
        ### Calculate monthly sums ###
        os.system("cdo " + fun + " " + files_to_merge1[file_ms] + " " + temp_dir_path + "/" + 
              variable + "_" + str(file_ms) + ".nc")

        ### Merge the monthly sum data ###

    #Rainfall is in mm/s, need to convert
    if variable=="pr":
        
        #temporary file
        temp_outfile = str(scratch_path + "/monthly_sums/" +"temp_precip" + 
                           bias_corr + "_" + model + "_" + scenario + ".nc")
        
        #Remove temporary file
        if os.path.isfile(temp_outfile): os.system("rm " + temp_outfile)

        
        os.system("cdo mergetime " + temp_dir_path + "/" + variable + "*.nc " + 
                  temp_outfile)

        #Multiply by the number of days and seconds per day
        os.system("cdo expr,'pr=pr*60*60*24' " + temp_outfile + " " +
        outfile)
        
        #Remove temporary file
        os.system("rm " + temp_outfile)
        
    else: 
        os.system("cdo mergetime " + temp_dir_path + "/" + variable + 
                  "*.nc " + outfile)

    
    
    
    
    
    
    
    
    
    
    
    


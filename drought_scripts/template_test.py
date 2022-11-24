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

lib_path  = "/g/data/w97/amu561/Steven_CABLE_runs/drought_scripts/functions"

##### ALTERTED #####
data_path = "/g/data/wj02/COMPLIANT"
#data_path = "/g/data/wj02/AWRA_OUTPUT"

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
perc=15
##########################################
##########################################
##########################################
##########################################
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
baseline=[1960,2005]

##########################
### FILE PREPROCESSING ###
##########################

### Define file locations ###
input_variable=['pr']
output_variable=['qtot','s0']

compliant=['QME','CCAM','MRNBC','ISIMIP2b']
#non_compliant=['RAW-GCM','NOBC-CCAM']

# Any output data from CCAM has a slightly differenet path
ccam_add=""
if bias_corr=="CCAM":
    ccam_add="CSIRO-CCAM-r3355-"

awra_add=""
if variable in output_variable:
    awra_add='AWRALv6-1-'

if variable in input_variable and bias_corr in compliant:
    data_path_var= data_path + "/HMINPUT/output/AUS-5/BoM"

if variable in output_variable and bias_corr in compliant:
    data_path_var= data_path + "/HMOUTPUT/output/AUS-5/BoM"

### Get historical and future simulations ###

if bias_corr=="CCAM":
    files_1_string=f'{data_path_var}/{awra_add}{model}/historical/r1i1p1/{ccam_add}r240x120-ISIMIP2b-AWAP/latest/day/{variable}/*.nc'
else:
    files_1_string=f'{data_path_var}/{awra_add}{model}/historical/r1i1p1/{ccam_add}r240x120-{bias_corr}-AWAP/latest/day/{variable}/*.nc'
files_to_merge1=glob.glob(files_1_string)

if bias_corr=="CCAM":
 files_2_string=f'{data_path_var}/{awra_add}{model}/{scenario}/r1i1p1/{ccam_add}r240x120-ISIMIP2b-AWAP/latest/day/{variable}/*.nc'
else:
    files_2_string=f'{data_path_var}/{awra_add}{model}/{scenario}/r1i1p1/{ccam_add}r240x120-{bias_corr}-AWAP/latest/day/{variable}/*.nc'
files_to_merge2=glob.glob(files_2_string)

files_to_merge1.extend(files_to_merge2)


### Create temporary directory ###
temp_dir_path = f"/scratch/w97/amu561/temp/{bias_corr}{model}{scenario}{variable}{str(scale)}"
os.system("mkdir -p " + temp_dir_path)

for file_ms in range(len(files_to_merge1)):
    
    ### Calculate monthly sums ###
    os.system("cdo monsum " + files_to_merge1[file_ms] + " " + temp_dir_path + "/" + 
              variable + "_" + str(file_ms) + ".nc")

### Merge the monthly sum data ###
os.system("cdo mergetime " + temp_dir_path + "/" + variable + "*.nc /scratch/w97/amu561/monthly_sums/"+
          bias_corr + "_" + model + "_" + scenario + "_" + variable + ".nc")

### Location of output file ###
files= str("/scratch/w97/amu561/outputs/monthly_sums/" + bias_corr + "_" + 
           model + "_" + scenario + "_" + variable + ".nc")

### Get ss data to add to s0 ###

if variable=="s0":

    ### Define file locations ###
    
    xtr1_files_to_merge=glob.glob(f'{data_path_var}/{awra_add}{model}/historical/r1i1p1/{ccam_add}r240x120-{bias_corr}-AWAP/latest/day/ss/*.nc')
    
    xtr2_files_to_merge=glob.glob(f'{data_path_var}/{awra_add}{model}/{scenario}/r1i1p1/{ccam_add}r240x120-{bias_corr}-AWAP/latest/day/ss/*.nc')

    xtr1_files_to_merge.extend(xtr2_files_to_merge)
    
    for file_ms in range(len(xtr1_files_to_merge)):

        ### Calculate monthly sums ###
        os.system("cdo monsum " + xtr1_files_to_merge[file_ms] + " " + temp_dir_path +
                  "/ss_" + str(file_ms) + ".nc")

    ### Merge monthly sum files ###
    os.system("cdo mergetime " + temp_dir_path + "/ss_*.nc /scratch/w97/amu561/monthly_sums/" +
              bias_corr + "_" + model + "_" + scenario + "_ss.nc")

    ### Location of output file ###
    xtr_files = str("/scratch/w97/amu561/monthly_sums/" + bias_corr + "_" + model + "_" + 
                    scenario + "_ss.nc")

### Delete temporary directory ###
#os.system("rm -r " + temp_dir_path) 

#################
### Load data ###
#################

#Model data
if variable == "s0":

    
    fh        = Dataset(files, mode='r')
    ds_ss     = Dataset(xtr_files, mode='r')
    s0_data   = fh.variables["s0"][:]
   
    ss_data   = ds_ss.variables["ss"][:]

    all_data  = s0_data.data + ss_data.data
    
    data      = all_data
    fh_time   = fh.variables["time"]
    mask      = s0_data.mask


else:
    
    fh       = Dataset(files, mode='r')
    all_data = fh.variables[variable][:] #[yr_ind]
    
    data     = all_data.data 
    mask     = all_data.mask
    fh_time  = fh.variables["time"]

#Get lon and lat (name varies by CMIP5 model)
try:
    lat = fh.variables['latitude'][:]
    lon = fh.variables['longitude'][:]
except:
    lat = fh.variables['lat'][:] #northing
    lon = fh.variables['lon'][:] #easting

#Get dataset dates
fh_dates = num2date(fh_time[:], fh_time.units, calendar=fh_time.calendar)
fh_years = np.array([y.year for y in fh_dates])

miss_val = -999
data[mask==True] = miss_val

control_ref = data

###################################################
### Create output file name and check if exists ###
###################################################
     
### Define output path ###
out_path = f'/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics/{scale}-month/{bias_corr}/{model}'

if not os.path.exists(out_path):    
    os.makedirs(out_path)
##########################################
# ##########################################
# ##########################################
# ##########################################       
#Create output file name
if variable=='s0':
    var_name='sm'
else:
    var_name=variable
out_file = str(out_path + "/drought_metrics_" + bias_corr + "_" + model + "_" + var_name + "_" + 
        scenario + "_" + str(scale) + ".nc")
##########################################
##########################################
##########################################
##########################################
##########################################


#############################
### Find baseline indices ###
#############################

#Get dates
ref_years = fh_years
if obs_ref:
    ref_years = obs_years
    

subset = range(np.where(ref_years == baseline[0])[0][0],
                np.where(ref_years == baseline[1])[0][-1] + 1) #Stupid python indexing
                
################################
### Initialise output arrays ###
################################

if return_all_tsteps:
    save_len = len(data)
else:
    save_len = int(len(data)*(perc/100)*2)

duration      = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
rel_intensity = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
intensity     = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
timing        = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    
tseries       = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    

if monthly:
    threshold    = np.zeros((12, len(lat), len(lon))) + miss_val # * np.nan
else:
    threshold    = np.zeros((len(lat), len(lon))) + miss_val # * np.nan

#########################
### Calculate metrics ###
#########################

#Loop through grid cells
for i in range(len(lat)):
            
    for j in range(len(lon)):
    
            #Calculate metrics if cell not missing
            if any(~mask[:,i,j]):  
            
                #Calculate metrics
                metric = drought_metrics(mod_vec=data[:,i,j], lib_path=lib_path, perc=perc, 
                                        monthly=monthly, obs_vec=control_ref[:,i,j],
                                        return_all_tsteps=return_all_tsteps, scale=scale,
                                        add_metrics=(['timing', 'rel_intensity', 'intensity', 'threshold']),
                                        subset=subset)
        
                ### Write metrics to variables ###
                duration[range(np.size(metric['duration'])),i,j]   = metric['duration']  #total drought duration (months)

                rel_intensity[range(np.size(metric['rel_intensity'])),i,j] = metric['rel_intensity'] #average magnitude
            
                intensity[range(np.size(metric['intensity'])),i,j] = metric['intensity'] #average intensity
    
                timing[range(np.size(metric['timing'])),i,j]       = metric['timing']    #drought timing (month index)

                tseries[range(np.size(metric['tseries'])),i,j]       = metric['tseries']    #drought timing (month index)

                if monthly:
                    threshold[:,i,j] = metric['threshold'][0:12]    #drought timing (month index)
                else:
                    threshold[i,j]   = metric['threshold']


##############################
### Write result to NetCDF ###
##############################
            
# Open a new netCDF file for writing
ncfile = Dataset(out_file,'w', format="NETCDF4_CLASSIC") 

# Create the output data
# Create the x, y and time dimensions
ncfile.createDimension('lat', lat.shape[0])
ncfile.createDimension('lon', lon.shape[0])
ncfile.createDimension('time', save_len)
    
if monthly:
    ncfile.createDimension('month', 12)

# Create dimension variables

longitude = ncfile.createVariable("lon",  'f8', ('lon',))
latitude  = ncfile.createVariable("lat",  'f8', ('lat',))
time      = ncfile.createVariable("time", 'i4', ('time',))

if monthly:
    month = ncfile.createVariable("month", 'i4', ('month',))

#Create data variables
data_dur  = ncfile.createVariable('duration', 'f8',('time','lat','lon'), fill_value=miss_val)
data_mag  = ncfile.createVariable('rel_intensity','f8',('time','lat','lon'), fill_value=miss_val)
data_int  = ncfile.createVariable('intensity','f8',('time','lat','lon'), fill_value=miss_val)
data_tim  = ncfile.createVariable('timing',   'i4',('time','lat','lon'), fill_value=miss_val)
data_ts   = ncfile.createVariable('tseries',   'i4',('time','lat','lon'), fill_value=miss_val)

#Create data variable for threshold
if monthly:
    data_thr = ncfile.createVariable('threshold', 'f8',('month','lat','lon'), fill_value=miss_val)
else:
    data_thr = ncfile.createVariable('threshold', 'f8',('lat','lon'), fill_value=miss_val)


#Set variable attributes
longitude.units = 'degrees_east'
latitude.units  = 'degrees_north'
time.units      = 'days since 1900-01-01'

time.calendar   = 'gregorian'

data_dur.long_name = 'drought event duration (no. months)'
data_mag.long_name = 'drought event relative intensity (%)'
data_int.long_name = 'drought event intensity (mm)'
data_tim.long_name = 'drought event timing (month index)'
data_thr.long_name = 'drought threshold (mm)'
data_ts.long_name  = 'original time series'

if monthly:
    month[:] = range(1,12+1)

# Write data to dimension variables
longitude[:]=lon
latitude[:] =lat

#If saving all time steps
if return_all_tsteps:
    time[:] = fh_time[:]
else:
    time[:] = range(1, save_len+1)
        
if monthly:
    month[:] = range(1,12+1)

#Write data to data variables
data_dur[:,:,:] = duration    
data_mag[:,:,:] = rel_intensity
data_int[:,:,:] = intensity
data_tim[:,:,:] = timing
data_ts[:,:,:]  = tseries

if monthly:    
    data_thr[:,:,:] = threshold
else:
    data_thr[:,:] = threshold

# Close the file
ncfile.close()


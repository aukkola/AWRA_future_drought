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
data_path = "/g/data/w97/amu561/Steven_CABLE_runs/"
#data_path = "/g/data/wj02/AWRA_OUTPUT"

scratch_path = '/scratch/w97/amu561/'


# Add lib_path to os directory
sys.path.append(os.path.abspath(lib_path))

#import all droughtmetric calcuation functions
from drought_metrics import *

########################
### DEFINE VARIABLES ###
########################


### Set Models ###
model="MIROC-MIROC5"

### Set variable ###
variable="qtot"

#CO2 directory
co2 = "CO2"

#Output file
out_file="test.nc" #str(sys.argv[4])


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
baseline=[1970,2005]

##########################
### FILE PREPROCESSING ###
##########################



### Get historical and future simulations ###

#historical
files_hist_string=str(data_path + '/CABLE_outputs/' + co2 +'_compressed/' + model + 
                   '/historical/r240x120-MRNBC-AWAP/outputs/' + '/*.nc')


#Files to merge                    
files_to_merge=glob.glob(files_hist_string)

#Future
files_fut_string=str(data_path + '/CABLE_outputs/' + co2 +'/' + model + 
                   '/rcp45/r240x120-MRNBC-AWAP/outputs/' + '/*.nc')


#Combine historical and future file names
files_to_merge_fut=glob.glob(files_fut_string)

files_to_merge.extend(files_to_merge_fut)


### Create temporary directory ###
temp_dir = str(scratch_path + "/monthly_sums_CABLE/")
os.system("mkdir -p " + temp_dir)

### Location of output file ###
temp_file= str(temp_dir + "/" + co2 + "_" + 
           model + "_" + variable + ".nc")



#If output file doesn't alreade exist, create it. Else skip 
if not os.path.isfile(temp_file):
    
    ds = xr.open_mfdataset(files_to_merge)

    if variable == "qtot" :
        
        total_q = ds.Qsb + ds.Qs

        #Multiple by the number of days in the month, and number of seconds per day
        #to convert to mm/month
        data = total_q * total_q.time.dt.daysinmonth * 86400.0


    elif variable== "sm":
        
        sm_all = ds.SoilMoist
        
        #Create weights by soil layer
        
        #Soil layer depths
        zse  = np.array([0.022, 0.058, 0.154, 0.409, 1.085, 2.872])
        weights_vec = zse / sum(zse)
        
        #Turn into xarray data array (otherwise xarray not happy)
        weights=xr.DataArray(weights_vec, dims=['soil'] )
        weights.name = "weights"
        
        #Take the weighted mean over soil layers to get one soil moisture value
        data = sm_all.weighted(weights).mean(["soil"])


    ### Write to file ###
        
    #Rename data variable
    data.name = variable
    
    #Write to file
    data.to_netcdf(temp_file, format='NETCDF4', 
                   encoding={variable:{
                             'shuffle':True,
                             'chunksizes':[12, 681, 40],
                             'zlib':True,
                             'complevel':5}
                             })


### Delete temporary directory ###
#os.system("rm -r " + temp_dir_path) 

#################
### Load data ###
#################

fh       = Dataset(temp_file, mode='r')
all_data = fh.variables[variable][:]

data     = all_data.data 
fh_time  = fh.variables["time"]
mask     = all_data.mask


#Get lon and lat (name varies by CMIP5 model)
try:
    lat = fh.variables['y'][:]
    lon = fh.variables['x'][:]
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
     
# ### Define output path ###
# out_path = f'/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics_CABLE/{scale}-month/{co2}/{model}/'
# 
# if not os.path.exists(out_path):    
#     os.makedirs(out_path)
# 
# #Create output file name
# out_file = str(out_path + "/drought_metrics_CABLE_"  + model + "_" + 
#                variable + "_rcp45_" + str(scale) + ".nc")
# 

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
    save_len = len(fh_dates)
else:
    save_len = int(len(data)*(perc/100)*2)

duration          = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
rel_intensity     = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
rel_intensity_mon = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan

#intensity     = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
timing        = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    
#tseries       = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    

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
                                        add_metrics=(['rel_intensity', 'threshold', 'rel_intensity_monthly',
                                        'timing']),
                                        subset=subset, miss_val=miss_val)
        
                ### Write metrics to variables ###
                duration[range(np.size(metric['duration'])),i,j]   = metric['duration']  #total drought duration (months)

                rel_intensity[range(np.size(metric['rel_intensity'])),i,j] = metric['rel_intensity'] #average magnitude

                rel_intensity_mon[range(np.size(metric['rel_intensity_monthly'])),i,j] = metric['rel_intensity_monthly'] #average magnitude
            
                #intensity[range(np.size(metric['intensity'])),i,j] = metric['intensity'] #average intensity
    
                timing[range(np.size(metric['timing'])),i,j]       = metric['timing']    #drought timing (month index)

                #tseries[range(np.size(metric['tseries'])),i,j]       = metric['tseries']    #drought timing (month index)

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
time      = ncfile.createVariable("time", 'f8', ('time',))

if monthly:
    month = ncfile.createVariable("month", 'i4', ('month',))

#Create data variables
data_dur  = ncfile.createVariable('duration', 'f8',('time','lat','lon'), fill_value=miss_val)
data_mag  = ncfile.createVariable('rel_intensity','f8',('time','lat','lon'), fill_value=miss_val)
data_rel  = ncfile.createVariable('rel_intensity_by_month','f8',('time','lat','lon'), fill_value=miss_val)

#data_int  = ncfile.createVariable('intensity','f8',('time','lat','lon'), fill_value=miss_val)
data_tim  = ncfile.createVariable('timing',   'i4',('time','lat','lon'), fill_value=miss_val)
#data_ts   = ncfile.createVariable('tseries',   'i4',('time','lat','lon'), fill_value=miss_val)

#Create data variable for threshold
if monthly:
    data_thr = ncfile.createVariable('threshold', 'f8',('month','lat','lon'), fill_value=miss_val)
else:
    data_thr = ncfile.createVariable('threshold', 'f8',('lat','lon'), fill_value=miss_val)


#Set variable attributes
longitude.units = 'degrees_east'
latitude.units  = 'degrees_north'
time.units      = fh_time.units

time.calendar   = fh_time.calendar

data_dur.long_name = 'drought event duration (no. months)'
data_mag.long_name = 'drought event relative intensity (%)'
data_rel.long_name = 'drought month relative intensity (%)'

#data_int.long_name = 'drought event intensity (mm)'
data_tim.long_name = 'drought event timing (binary drought/non-drought index)'
data_thr.long_name = 'drought threshold (mm)'
#data_ts.long_name  = 'original time series'

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
data_rel[:,:,:] = rel_intensity_mon

#data_int[:,:,:] = intensity
data_tim[:,:,:] = timing
#data_ts[:,:,:]  = tseries

if monthly:    
    data_thr[:,:,:] = threshold
else:
    data_thr[:,:] = threshold

# Close the file
ncfile.close()


#Finally compress file
#os.system("nccompress -o " + out_file) 



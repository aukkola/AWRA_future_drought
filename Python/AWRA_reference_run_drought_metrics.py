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

#Calculates droughts metrics for AWRA historical reference runs

### Set variable ###
variable=['qtot', 'sm']

### Set drought metric conditions ###
return_all_tsteps=True

### Set percentile for drought threshold ###
perc=15

### Set scale for month aggregation ###
scale=3

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

#Setting this in case have duplicate time steps
os.system("export SKIP_SAME_TIME=1")

#Make temporary folder
temp_dir_path=str(scratch_path + "/monthly_sums_AWRA_ref/")
os.system("mkdir -p " + temp_dir_path)


#Loop through variables
for v in range(len(variable)):

    #Runoff
    if variable[v] == "qtot":
        
        files = str("/g/data/fj8/BoM/AWRA/DATA/AWRA_REF_FORECAST_HYDROPROJ/" +
                              "awral_orv6qes-viney-icc.2018.1.163/sim/" + variable[v] + "_*.nc")
                  
        
        ### Merge and then calculate monthly sum data ###
        
        temp_qtot_file = str(temp_dir_path + "/qtot_temp.nc")
        
        os.system("cdo mergetime " + files + " " + temp_qtot_file)

        ### Calculate monthly sums ###
        os.system("cdo monsum -selyear,1960/2020 " + temp_qtot_file + " " + temp_dir_path +
                  "/qtot.nc")

        os.system("rm "+ temp_qtot_file)
    
                  
    #Soil moisture              
    elif variable[v] == "sm":
        
        files_s0 = str("/g/data/fj8/BoM/AWRA/DATA/AWRA_REF_FORECAST_HYDROPROJ/" +
                                 "awral_orv6qes-viney-icc.2018.1.163/sim/s0_*.nc")
                                 
        files_ss = str("/g/data/fj8/BoM/AWRA/DATA/AWRA_REF_FORECAST_HYDROPROJ/" +
                                 "awral_orv6qes-viney-icc.2018.1.163/sim/ss_*.nc")
        
        
        ### Merge and then calculate monthly sum data ###
        
        #s0 data
        
        temp_s0_file = str(temp_dir_path + "/s0_temp.nc")
        
        os.system("cdo mergetime " + files_s0 + " " + temp_s0_file)

        os.system("cdo monmean -selyear,1960/2020 " + temp_s0_file + " " + temp_dir_path +
                  "/s0.nc")

        os.system("rm "+ temp_s0_file)
    
        #ss data
        temp_ss_file = str(temp_dir_path + "/ss_temp.nc")
        
        os.system("cdo mergetime " + files_ss + " " + temp_ss_file)

        os.system("cdo monmean -selyear,1960/2020 " + temp_ss_file + " " + temp_dir_path +
                  "/ss.nc")

        os.system("rm "+ temp_ss_file)


    #################
    ### Load data ###
    #################

    #Model data
    if variable[v] == "sm":

        fh        = Dataset(str(temp_dir_path +"/s0.nc"), mode='r')
        ds_ss     = Dataset(str(temp_dir_path +"/ss.nc"), mode='r')
        s0_data   = fh.variables["s0"][:]
       
        ss_data   = ds_ss.variables["ss"][:]

        #Sum ss and s0
        data  = s0_data.data + ss_data.data
        
        fh_time   = fh.variables["time"]
        mask      = s0_data.mask


    elif variable[v] == "qtot":
        
        fh       = Dataset(str(temp_dir_path +"/qtot.nc"), mode='r')
        all_data = fh.variables[variable[v]][:] #[yr_ind]
        
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
    out_path = f'/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics_AWRA_ref/{scale}-month/'

    if not os.path.exists(out_path):    
        os.makedirs(out_path)
    ##########################################
    # ##########################################
    # ##########################################
    # ##########################################       
    #Create output file name
    var_name=variable[v]

    out_file = str(out_path + "/drought_metrics_AWRA_ref_" + "_scale_" + str(scale) + "_1960_2020.nc")
                   
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
    time      = ncfile.createVariable("time", 'i4', ('time',))

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
    time.units      = 'days since 1960-01-01'

    time.calendar   = 'gregorian'

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


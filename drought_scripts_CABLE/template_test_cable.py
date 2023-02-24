












import xarray as xr
import matplotlib.pyplot as plt
import sys 
import os

#Need to first fix CABLE leap year calendars
#If this code fails, make sure to run the Bash/fix_CABLE_calendar.sh script first

lib_path  = "/g/data/w97/amu561/Steven_CABLE_runs/scripts/drought_scripts/functions"
# Add lib_path to os directory
sys.path.append(os.path.abspath(lib_path))

#import all droughtmetric calcuation functions
from drought_metrics import *


path="/g/data/w97/amu561/Steven_CABLE_runs/CABLE_outputs/"
ds = xr.open_mfdataset(path + "CO2/MIROC-MIROC5/rcp45/r240x120-MRNBC-AWAP/outputs/*.nc")


total_q = ds.Qsb + ds.Qs


#Multiple by the number of days in the month, and number of seconds per day
#to convert to mm/month
total_q = total_q * total_q.time.dt.daysinmonth * 86400.0

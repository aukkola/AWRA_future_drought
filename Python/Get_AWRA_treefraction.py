
import h5py
from netCDF4 import Dataset,num2date


#Read AWRA parameter file
f = h5py.File('/g/data/er4/AWRACMS/awrams_cm_v6.1/models/awrams/models/awral/data/spatial_parameters.h5', 'r')


#Parameter fields
params = f['parameters']

#Get tree fraction field
ftree = params['f_tree']

#Dimension fields
dims = f['dimensions']

#Get longitude and latitude
lon = dims['longitude']
lat = dims['latitude']


### Write NetCDF file ###

outfile = '/g/data/w35/amu561/Steven_CABLE_runs/AWRA_fields/AWRA_tree_fraction.nc'            

# Open a new netCDF file for writing
ncfile = Dataset(outfile,'w', format="NETCDF4_CLASSIC") 

# Create the output data
# Create the x, y and time dimensions
ncfile.createDimension('lat', lat.shape[0])
ncfile.createDimension('lon', lon.shape[0])
    

# Create dimension variables
longitude = ncfile.createVariable("lon",  'f8', ('lon',))
latitude  = ncfile.createVariable("lat",  'f8', ('lat',))

#Data variable
data  = ncfile.createVariable('tree_fraction', 'f8',('lat','lon'), fill_value=-9999)

#Set variable attributes
longitude.units = 'degrees_east'
latitude.units  = 'degrees_north'

# Write data to dimension variables
longitude[:]=lon[:]
latitude[:] =lat[:]

#Write data to data variables
data[:,:] = ftree[:,:]

                
# Close the file
ncfile.close()



 









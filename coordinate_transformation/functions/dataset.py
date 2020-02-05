from pathlib import Path
from netCDF4 import Dataset
import numpy as np 
import datetime as dt

# Create .nc file
f = Dataset('topography_coord.nc','w', format='NETCDF4')
f.description = 'Togography data in sph, Cart, and cyl coordinates'

# Create dimensions
f.createDimension('colat', len(colat_Prt_cnt))
f.createDimension('lon', len(lon_Prt))

# Create variables, 'f4' for single precision floats, i.e. 32bit
radial_distance = f.createVariable('radial_distance', 'f4', ('colat', 'lon'))
colatitude = f.createVariable('colatitude', 'f4', 'colat')
longitude = f.createVariable('longitude', 'f4', 'lon')
x_value = f.createVariable('x_value', 'f4', ('colat', 'lon'))
y_value = f.createVariable('y_value', 'f4', ('colat', 'lon'))
z_value = f.createVariable('z_value', 'f4', ('colat', 'lon'))
cyl_radial_distance = f.createVariable('cyl_radial_distance', 'f4', ('colat', 'lon'))
azimuth = f.createVariable('azimuth', 'f4', 'lon')
height = f.createVariable('height', 'f4', ('colat', 'lon'))

radial_distance [:] = r_cnt_bathy
colatitude [:] = colat_Prt_cnt
longitude [:] = lon_Prt
x_value [:] = x_Prt
y_value [:] = y_Prt
z_value [:] = z_Prt
cyl_radial_distance [:] = s_Prt
azimuth [:] = phi_Prt
height [:] = l_Prt

# Add attributes to the file
today = dt.datetime.now()
f.history = "Created " + today.strftime("%d/%m/%y")
#Add local attributes to variable instances
radial_distance.units = 'km'
colatitude.units = 'degrees north'
longitude.units = 'degrees east'
x_value.units = 'km'
y_value.units = 'km'
z_value.units = 'km'
cyl_radial_distance.units = 'km'
azimuth.units = 'degrees east'
height.units = 'km'

f.close()

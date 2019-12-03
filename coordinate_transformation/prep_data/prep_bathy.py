import numpy as np
from math import pi, sin, cos, sqrt
from pathlib import Path
import netCDF4 as nc4
import numpy.ma as ma
import datetime as dt

from transformation_functions.get_domain import find_nearest, truncate_domain
from transformation_functions.get_spherical import wgs84, geograph_to_geocent, radius_cnt
from transformation_functions.get_cartesian import sph_to_cartesian
from transformation_functions.get_cylindrical import sph_to_cylindrical
from transformation_functions.get_rotation import rotation_matrix, rotate_N_pole

data_folder = Path('/Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/Code/prep_data/GEBCO_2019/')
file2open = data_folder / 'GEBCO_2019.nc' #file with location
nc_GEBCO = nc4.Dataset(file2open, 'r')
raw_lat = nc_GEBCO.variables['lat'][:] # in degrees N
raw_lon = nc_GEBCO.variables['lon'] [:]# in degrees E 
raw_elevation = nc_GEBCO.variables['elevation'] [:] 
# in m as height above reference ellipsoid

# Unmask arrays
raw_lat = np.ma.getdata(raw_lat)
raw_lon = np.ma.getdata(raw_lon)
raw_elevation = np.ma.getdata(raw_elevation)

#%% Truncate, transform, rotate
# Define domain
# Find indices for 35.5-41.4N & -22 - -14.5E 
lat_max = 41.4 
lat_min = 35.5
lon_max = -14.5
lon_min = -22
bounds = [lat_max, lat_min, lon_max, lon_min]

# Truncate domain
(lat_Prt, lon_Prt, bathy_Prt) = truncate_domain(raw_lat, raw_lon, raw_elevation, bounds)

# Transform latitudes from geographic to geocentric
lat_Prt_cnt = geograph_to_geocent(lat_Prt)

# Get spherical coordinates
len_lon = len(lon_Prt)
# Calculate radius of reference surface at geocentric latitudes
r_cnt = radius_cnt(lat_Prt_cnt) 
r_cnt = np.array([r_cnt,]*len_lon).conj().transpose()
# Calculate the radius for each bathymetry data point, in km 
r_cnt_bathy = (r_cnt + bathy_Prt)/1000 
# Calculate the colatitude to define a spherical coordinate system
colat_Prt_cnt = 90 - lat_Prt_cnt

# Get Cartesian
(x_Prt,y_Prt,z_Prt) = sph_to_cartesian(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

# Get cylindrical
(s_Prt,phi_Prt,l_Prt) = sph_to_cylindrical(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

# Rotate Cartesian & back-transform to spherical?
#Find source lat and lon (centre of domain), calculate length of degrees
#to do it accurately but at the moment just going to take the average of 
#lat and lon as the centre
av_lat = np.mean(lat_Prt) # in geographic
av_lon = np.mean(lon_Prt)
(x_rot, y_rot, z_rot) = rotate_N_pole(av_lat, av_lon, x_Prt, y_Prt, z_Prt)

#%% Save data as netCDF files

### Get .nc file for different coordinate systems
# Create .nc file
f = nc4.Dataset('topography_coord.nc','w', format='NETCDF4')
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
# Assign values to variables
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
# Close file
f.close()


### Get .nc file for rotated file
g = nc4.Dataset('bathymetry_N_pole.nc','w', format='NETCDF4')
g.description = 'Togography data Cartesian coordinates, rotated to be centred about N pole'
# Create dimensions
g.createDimension('colat', len(colat_Prt_cnt))
g.createDimension('lon', len(lon_Prt))
# Create variables
x_rot_value = g.createVariable('x_rot_value', 'f4', ('colat', 'lon'))
y_rot_value = g.createVariable('y_rot_value', 'f4', ('colat', 'lon'))
z_rot_value = g.createVariable('z_rot_value', 'f4', ('colat', 'lon'))
# Assign values to variables
x_rot_value [:] = x_rot
y_rot_value [:] = y_rot
z_rot_value [:] = z_rot
# Add attributes to the file
today = dt.datetime.now()
g.history = "Created " + today.strftime("%d/%m/%y")
#Add local attributes to variable instances
x_value.units = 'km'
y_value.units = 'km'
z_value.units = 'km'
# Close file
g.close()

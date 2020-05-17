# %%
import numpy as np
from math import pi, sin, cos, sqrt
from pathlib import Path
import numpy.ma as ma
import datetime as dt

from coordinate_transformation.functions.domain import find_nearest, truncate_domain
from coordinate_transformation.functions.get_rotation \
    import *
from coordinate_transformation.functions.domain \
    import get_variable

# %% read data
data_folder = Path('coordinate_transformation/raw_data/crust1.0/')
file2open = data_folder / 'depthtomoho.xyz' #file with location
# lon (deg E), lat (deg N), depth (km, negative down)
f = open(file2open, 'r')
contents = np.loadtxt(f, usecols=[0,1,2])
lon = contents[:,0]
lat = contents[:,1]
moho = contents[:,2]



# in m as height above or below?? reference ellipsoid

# # Unmask arrays
# raw_lat = np.ma.getdata(raw_lat)
# raw_lon = np.ma.getdata(raw_lon)
# raw_elevation = np.ma.getdata(raw_elevation)

# %% Truncate, transform, rotate
# # Define domain
# # Find indices for 35.5-41.4N & -22 - -14.5E 
# lat_max = 41.4 
# lat_min = 35.5
# lon_max = -14.5
# lon_min = -22
# bounds = [lat_max, lat_min, lon_max, lon_min]

# # Truncate domain
# (lat_Prt, lon_Prt, bathy_Prt) = truncate_domain(raw_lat, raw_lon, raw_elevation, bounds)

# # Transform latitudes from geographic to geocentric
# lat_Prt_cnt = geograph_to_geocent(lat_Prt)

# # Get spherical coordinates
# len_lon = len(lon_Prt)
# # Calculate radius of reference surface at geocentric latitudes
# r_cnt = radius_cnt(lat_Prt_cnt) 
# r_cnt = np.array([r_cnt,]*len_lon).conj().transpose()
# # Calculate the radius for each bathymetry data point, in km 
# r_cnt_bathy = (r_cnt + bathy_Prt)/1000 
# # Calculate the colatitude to define a spherical coordinate system
# colat_Prt_cnt = 90 - lat_Prt_cnt

# # Get Cartesian
# (x_Prt,y_Prt,z_Prt) = sph_to_cartesian(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

# # Get cylindrical
# (s_Prt,phi_Prt,l_Prt) = sph_to_cylindrical(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

# # Rotate Cartesian & back-transform to spherical?
# #Find source lat and lon (centre of domain), calculate length of degrees
# #to do it accurately but at the moment just going to take the average of 
# #lat and lon as the centre
# av_lat = np.mean(lat_Prt) # in geographic
# av_lon = np.mean(lon_Prt)
# (x_rot, y_rot, z_rot) = rotate_N_pole(av_lat, av_lon, x_Prt, y_Prt, z_Prt)

# %% Save data as netCDF files
# # Get .nc file for different coordinate systems
# # Create .nc file
# f = Dataset('topography_coord.nc','w', format='NETCDF4')
# f.description = 'Togography data in sph, Cart, and cyl coordinates'

# # Create dimensions
# f.createDimension('colat', len(colat_Prt_cnt))
# f.createDimension('lon', len(lon_Prt))

# # Create variables, 'f4' for single precision floats, i.e. 32bit
# radial_distance = f.createVariable('radial_distance', 'f4', ('colat', 'lon'))
# colatitude = f.createVariable('colatitude', 'f4', 'colat')
# longitude = f.createVariable('longitude', 'f4', 'lon')
# x_value = f.createVariable('x_value', 'f4', ('colat', 'lon'))
# y_value = f.createVariable('y_value', 'f4', ('colat', 'lon'))
# z_value = f.createVariable('z_value', 'f4', ('colat', 'lon'))
# cyl_radial_distance = f.createVariable('cyl_radial_distance', 'f4', ('colat', 'lon'))
# azimuth = f.createVariable('azimuth', 'f4', 'lon')
# height = f.createVariable('height', 'f4', ('colat', 'lon'))

# radial_distance [:] = r_cnt_bathy
# colatitude [:] = colat_Prt_cnt
# longitude [:] = lon_Prt
# x_value [:] = x_Prt
# y_value [:] = y_Prt
# z_value [:] = z_Prt
# cyl_radial_distance [:] = s_Prt
# azimuth [:] = phi_Prt
# height [:] = l_Prt

# # Add attributes to the file
# today = dt.datetime.now()
# f.history = "Created " + today.strftime("%d/%m/%y")
# #Add local attributes to variable instances
# radial_distance.units = 'km'
# colatitude.units = 'degrees north'
# longitude.units = 'degrees east'
# x_value.units = 'km'
# y_value.units = 'km'
# z_value.units = 'km'
# cyl_radial_distance.units = 'km'
# azimuth.units = 'degrees east'
# height.units = 'km'

# f.close()


# # Get .nc file for rotated file

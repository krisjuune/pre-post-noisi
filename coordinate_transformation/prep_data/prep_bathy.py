# %% Get raw data
import numpy as np
from math import pi, sin, cos, sqrt
from pathlib import Path
import netCDF4 as nc4
import numpy.ma as ma
import datetime as dt

from coordinate_transformation.functions.get_domain import \
    find_nearest, truncate_domain
from coordinate_transformation.functions.get_spherical import \
    wgs84, geographic_to_geocentric, radius_cnt
from coordinate_transformation.functions.get_cartesian import \
    sph_to_cartesian
from coordinate_transformation.functions.get_cylindrical import \
    sph_to_cylindrical
from coordinate_transformation.functions.get_rotation import \
    rotation_matrix, rotate_N_pole

data_folder = Path('coordinate_transformation/raw_data/GEBCO_2019')
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

# %% Truncate, transform, rotate
from coordinate_transformation.functions.get_rotation \
    import get_cartesian_distance
# Define domain
# Find indices for 35.5...39.5N, -14..-19E 
# Made slightly bigger below to cover for sure the whole area with topography
lat_max = 40
lat_min = 35
lon_max = -13.5
lon_min = -19.5
bounds = [lat_max, lat_min, lon_max, lon_min]

# Truncate domain
(lat_Prt, lon_Prt, bathy_Prt) = truncate_domain(raw_lat, raw_lon, raw_elevation, bounds)

# Transform latitudes from geographic to geocentric
lat_Prt_cnt = geographic_to_geocentric(lat_Prt)

# Get lat and lon as distances from the centre of domain (N pole)
(x_N, y_N) = get_cartesian_distance(lon_Prt, lat_Prt_cnt)

# Express bathymetry as depth from reference level (4.72km) in km
# positive downwards (i.e. deeper)

rel_bathymetry = 4.72*(-1) - bathy_Prt/1000

# %% Save data as netCDF files

### Get .nc file for the x, y, and depth variables
# Create .nc file
f = nc4.Dataset('topography_coord.nc','w', format='NETCDF4')
f.description = 'Togography data in Cartesian coordinates'
# Create dimensions
f.createDimension('x', len(x_N))
f.createDimension('y', len(y_N))
# Create variables, 'f4' for single precision floats, i.e. 32bit
x = f.createVariable('x', 'f4', 'x')
y = f.createVariable('y', 'f4', 'y')
z = f.createVariable('bathymetry', 'f4', ('x', 'y'))
# Assign values to variables
x [:] = x_N
y [:] = y_N
z [:] = rel_bathymetry
# Add attributes to the file
today = dt.datetime.now()
f.history = "Created " + today.strftime("%d/%m/%y")
#Add local attributes to variable instances
x.units = 'km'
y.units = 'km'
z.units = 'km'
# Close file
f.close()


# %%

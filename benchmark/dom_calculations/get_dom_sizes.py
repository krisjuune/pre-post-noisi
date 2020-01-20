# %% Load variables in WGS84

import pickle
import numpy as np
import numpy.ma as ma

with open('coordinate_transformation/variables/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)
with open('coordinate_transformation/variables/colat_Prt_cnt', 'rb') as f:
    colat_Prt_cnt = pickle.load(f)
colat_Prt_cnt = np.ma.getdata(colat_Prt_cnt)
lat_Prt_cnt = 90 - colat_Prt_cnt
with open('coordinate_transformation/variables/bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)

# %% Truncate domain for mesh validation 
from coordinate_transformation.transformation_functions.get_domain import truncate_domain, find_nearest

# Define domain, pretty arbitrarily chosen 
# To be able to calculate sizes lat_min and max
# must be integers.
lat_max = 40
lat_min = 36
lon_max = -15
lon_min = -20.5
bounds = [lat_max, lat_min, lon_max, lon_min]

# Truncate domain
(lat_dom, lon_dom, r_dom) = truncate_domain(lat_Prt_cnt, lon_Prt, bathy_Prt, bounds)

# %% Calcualte size of domain in km for simulation nr 3
import numpy as np
from math import pi
from coordinate_transformation.transformation_functions.get_spherical \
    import wgs84
from benchmark.processing.processing_functions \
    import len_deg_lat, len_deg_lon

lat_int = np.arange(lat_min,lat_max, dtype = int)
dlon = len_deg_lon(lat_int)
dlat = len_deg_lat(lat_int)

# in km
range_lon_min = (lon_max - lon_min)*np.amin(dlon)/1000 
range_lat = np.sum(dlat)/1000

print('Geographic radius is', range_lon_min/2, \
    'along lon at max lat, in km', 'for bounds', bounds)
print('Geographic radius is', range_lat/2, \
    'along lat, in km', 'for bounds', bounds)

#%% Get Cartesian equivalent for simulation nr 1
import numpy as np
from math import pi, sin, cos, sqrt
from transformation_functions.get_spherical \
    import wgs84, geograph_to_geocent
# Assume spherical Earth, and take the 
# Cartesian radius to be the far cathetus of 
# the triangle forming between the radius at 
# the centre of the domain, the drawn tangent, 
# and the minimum (or max, doesn't matter) 
# latitude. 

# Calculate radius at the centre of the domain
lat_mid = (lat_max + lat_min)/2
lat_mid = geograph_to_geocent(lat_mid)
lat_min_cnt = geograph_to_geocent(lat_min)
a = wgs84()[0]
b = wgs84()[1]
# Radius at the centre of domain, i.e. at 38N
r38 = np.sqrt((a**2*(np.cos(np.deg2rad(lat_mid))**2)) + \
    (b**2*(np.sin(np.deg2rad(lat_mid))**2)))
# Calculate Cartesian equivalent range, in km
range_Cart = r38*np.sin(np.deg2rad(lat_mid - \
    lat_min_cnt))/1000
print('Cartesian eqv radius is ', range_Cart, \
    'km for bounds', bounds)

# %% Get spherical equivalent for simulation nr 2
# Assume spherical Earth with the radius defined 
# by the radius at the centre of the domain

range_Sph_lat = (lat_max-lat_min)*len_deg_lat(lat_mid)/1000
range_Sph_lon = (lon_max-lon_min)*len_deg_lon(lat_mid)/1000
print('Spherical eqv radius is ', \
    range_Sph_lon/2, \
    'along lon at mid lat, in km', \
    'for bounds', bounds)
print('Spherical eqv radius is ', \
    range_Sph_lat/2, \
    'along lat at mid lat, in km', \
    'for bounds', bounds)

# %%

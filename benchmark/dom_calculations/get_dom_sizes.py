# %% Load variables in WGS84

import pickle
import numpy as np
import numpy.ma as ma

with open('coordinate_transformation/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)
with open('coordinate_transformation/colat_Prt_cnt', 'rb') as f:
    colat_Prt_cnt = pickle.load(f)
colat_Prt_cnt = np.ma.getdata(colat_Prt_cnt)
lat_Prt_cnt = 90 - colat_Prt_cnt
with open('coordinate_transformation/bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)

# %% Truncate domain for mesh validation 
from transformation_functions.get_domain import truncate_domain, find_nearest

# Define domain, pretty arbitrarily chosen 
lat_max = 39
lat_min = 37
lon_max = -16.5
lon_min = -19
bounds = [lat_max, lat_min, lon_max, lon_min]

# Truncate domain
(lat_dom, lon_dom, r_dom) = truncate_domain(lat_Prt_cnt, lon_Prt, bathy_Prt, bounds)

# %% Calcualte size of domain in km for simulation nr 3
import numpy as np
from math import pi
# from transformation_functions.get_spherical import wgs84

lat_max = 39
lat_min = 37
lon_max = -16.5
lon_min = -19.5
bounds = [lat_max, lat_min, lon_max, lon_min]

def wgs84(): #WGSS84 coordinate system with Greenwich as lon = 0
    # set semi-major axis of the oblate spheroid Earth, in m
    a = 6378137.0
    # set semi-minor axis of the oblate spheroid Earth, in m
    b = 6356752.314245
    # calculate inverse flattening f
    f = a/(a-b)
    # calculate squared eccentricity e
    e_2 = (a**2-b**2)/a**2
    return(a,b,e_2,f)

# Calculate the length of one degree of lat and lon as a function of lat
def len_deg_lon(lat):
    e_2 = wgs84()[2]
    a = wgs84() [0]
    # This is the length of one degree of longitude 
    # approx. after WGS84, at latitude lat
    # in m
    lat = pi/180*lat
    dlon = (pi*a*np.cos(lat))/180*np.sqrt((1-e_2*np.sin(lat)**2))
    return np.round(dlon,5)

def len_deg_lat(lat):
    # This is the length of one degree of latitude 
    # approx. after WGS84, between lat-0.5deg and lat+0.5 deg
    # in m
    lat = pi/180*lat
    dlat = 111132.954 - 559.822 * np.cos(2*lat) + 1.175*np.cos(4*lat)
    return np.round(dlat,5)

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
# from transformation_functions.get_spherical import wgs84, \
    # geograph_to_geocent
# Assume spherical Earth, and take the Cartesian radius to be 
# the far cathetus of the triangle forming between the radius
# at the centre of the domain, the drawn tangent, and the 
# minimum (or max, doesn't matter) latitude. 

def geograph_to_geocent(theta):
    e_2 = wgs84()[2] # returns the 3rd item from function ouputs 
    theta = np.rad2deg(np.arctan((1 - e_2) * np.tan(np.deg2rad(theta))))
    return theta

# Calculate radius at the centre of the domain
lat_mid = (lat_max + lat_min)/2
lat_mid = geograph_to_geocent(lat_mid)
lat_min_cnt = geograph_to_geocent(lat_min)
a = wgs84()[0]
b = wgs84()[1]
r38 = np.sqrt((a**2*(np.cos(np.deg2rad(lat_mid))**2)) + \
    (b**2*(np.sin(np.deg2rad(lat_mid))**2)))
# Calculate Cartesian equivalent range, in km
range_Cart = r38*np.sin(np.deg2rad(lat_mid - lat_min_cnt))\
    /1000
print('Cartesian eqv radius is ', range_Cart, 'km for bounds', bounds)

# %% Get spherical equivalent for simulation nr 2
# Assume spherical Earth with the radius defined by the 
# radius at the centre of the domain

range_Sph_lat = (lat_max-lat_min)*len_deg_lat(lat_mid)/1000
range_Sph_lon = (lon_max-lon_min)*len_deg_lon(lat_mid)/1000
print('Spherical eqv radius is ', range_Sph_lon/2, 'along lon at mid lat, in km', \
    'for bounds', bounds)
print('Spherical eqv radius is ', range_Sph_lat/2, 'along lat at mid lat, in km', \
    'for bounds', bounds)

# %%

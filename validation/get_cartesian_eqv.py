#%% Load variables in WGS84

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

#%% Truncate domain for mesh validation 
from transformation_functions.get_domain import truncate_domain, find_nearest

# Define domain, pretty arbitrarily chosen 
lat_max = 39
lat_min = 37
lon_max = -16.5
lon_min = -19
bounds = [lat_max, lat_min, lon_max, lon_min]

# Truncate domain
(lat_dom, lon_dom, r_dom) = truncate_domain(lat_Prt_cnt, lon_Prt, bathy_Prt, bounds)

#%% Calcualte size of domain in km 
import numpy as np
from math import pi
from transformation_functions.get_spherical import wgs84

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

print(range_lon_min, 'along lon at max lat, in km', 'for bounds', bounds)
print(range_lat, 'along lat, in km', 'for bounds', bounds)

#%% Get average bathymetry and depth to Moho 
from pathlib import Path
# Could weight the average calculations by the len of lat and lon at each point 
# But do it here in a rubbish manner by just taking the unweighted average
print(np.mean(bathy_Prt), 'is the average bathymetry in m')

data_folder = Path('coordinate_transformation/raw_data/crust1.0/')
moho_file = data_folder / 'depthtomoho.xyz'
raw_moho = open(moho_file, 'r')
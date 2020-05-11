#%% import ... 

import numpy as np 
import numpy.ma as ma
import pandas as pn 
import pickle
import plotly.graph_objects as go 
from math import pi, cos, sin, sqrt
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from transformation_functions.get_spherical import wgs84

with open('coordinate_transformation/r_cnt_bathy', 'rb') as f:
    r_cnt_bathy = pickle.load(f)
r_cnt_bathy = np.ma.getdata(r_cnt_bathy)
with open('coordinate_transformation/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)
with open('coordinate_transformation/colat_Prt_cnt', 'rb') as f:
    colat_Prt_cnt = pickle.load(f)
colat_Prt_cnt = np.ma.getdata(colat_Prt_cnt)


# %% Plot spherical 
# Variables in spherical: r_cnt_bathy, colat_Prt_cnt, lon_Prt

# Calculate length of degree of lat and lon, by Laura Ermert
def len_deg_lon(lat):
    e_2 = wgs84()[2]
    a = wgs84() [0]
    # This is the length of one degree of longitude approximated
    # after WGS84 at latitude lat
    # in m
    lat = pi/180*lat
    dlon = (pi*a*cos(lat))/180*sqrt((1-e_2*sin(lat)**2)) 
    # Error only length 1 arrays can be converted to Python 
    # scalars
    return round(dlon,5)

def len_deg_lat(lat):
    # This is the length of one degree of latitude approximated
    # after WGS84, between lat-0.5deg and lat+0.5 deg
    # in m
    lat = pi/180*lat
    dlat = 111132.954 - 559.822 * cos(2*lat) + 1.175*cos(4*lat)
    return round(dlat,5)

# length_lon = len_deg_lon(lat_Prt_cnt)
# How to make a mesh?

# Create 3D surface plot
fig = go.Figure(data=[go.Surface(z=r_cnt_bathy)])
# in a scrappy manner (don't account for curvature, assuming 
# all degrees are same length)
# x - lon, y - lat, z - radial distance
fig.update_layout(autosize=True)
fig.show()
# checked with the scrappy Matlab figure pre transforming into 
# spherical coordinate system and it's the same, cry

#%% Plot spherical with geopandas
import pandas as pd
import geopandas as gpd

# %% Plot Cartesian pre rotation
# Variables in Cartesian pre-rotation: x_Prt, y_Prt, z_Prt

# Create 3D matplotlib plot
fig = plt.figure()
ha = fig.add_subplot(111, projection='3d')
ha.plot_surface(x_Prt, y_Prt, z_Prt)
plt.show()
# Doesn't work, how to get the proper x and y and depth in z?

# %% Plot Cartesian post rotation

# Doesn't work with matplotlib
# fig = plt.figure()
# ha = fig.add_subplot(111, projection='3d')
# ha.plot_surface(x_rot, y_rot, z_rot)
# plt.show()

fig = go.Figure(data=[go.Surface(z=z_rot)])
fig.update_layout(autosize=True)
fig.show()

# Also doesn't work with the surface plot but a profile of 
# topography is somehow there in the data, shows no change 
# with y and a weird curvature in x
# Need to also plot the x and y data but how??
# %%
# The Cartesian plots are clearly wrong, how to plot three
# arrays, each in 2D (1416x1800) in 3D, basically as a 
# surface plot???
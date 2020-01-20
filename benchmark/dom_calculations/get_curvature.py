# %% Load variables
import numpy as np 
import numpy.ma as ma
import pickle 
from pathlib import Path
from math import pi
from coordinate_transformation.transformation_functions.get_domain \
    import find_nearest, truncate_domain

with open('coordinate_transformation/variables/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)

with open('coordinate_transformation/variables/lat_Prt', 'rb') as f:
    lat_Prt = pickle.load(f)
lat_Prt = np.ma.getdata(lat_Prt)

lat_min = 36
lat_max = 40
lon_min = -20.5
lon_max = -15
bounds = [lat_max, lat_min, lon_max, lon_min]

(lat_Prt, lon_Prt) = truncate_domain(lat_Prt, \
    lon_Prt, bounds)

# %% Get arrays with distances to the curve

def get_curvature(lat, lon, radius = 6371, \
    theta = 38, phi = -17.75):
    """
    Function to calculate the curvature 
    relative to a flat surface at the given 
    radius assuming a sphere with the given 
    radius. Inputs include arrays of latitude, 
    longitude, and radius. Function returns 
    array of depths relative to this flat 
    surface with the dimensions of lon, lat.
    Units in degrees for angles, km for 
    distances.  
    """
    # preallocate output array
    curvature = np.zeros((len(lon), \
        len(lat)), float) 
    
    # convert to radians
    lon = pi/180*lon
    lat = pi/180*lat
    phi = pi/180*phi
    theta = pi/180*theta

    # loop over the lats and lons
    for i in range(len(lon)):
        for j in range(len(lat)):
            # find angle between point i,j and 
            # centre
            a = radius*np.sin(lon[i] - phi)
            b = radius*np.sin(lat[j] - theta)
            c = np.sqrt(np.square(a) + \
                np.square(b))
            # arcsin(x), x has to be [-1,1]
            alpha = np.arcsin(c/radius)
            # calculate depth to curve from flat 
            # surface
            y = radius/np.cos(alpha) - radius
            x = y*np.cos(alpha)
            curvature [i,j] = x 
    
    return(curvature)

spherical_curvature = get_curvature(lat_Prt, \
    lon_Prt)

# %% Get curvature for ellipsoid

def get_curvature_wgs84(lat, lon, radius):
    """
    Function to calculate the curvature 
    relative to a flat surface at the given 
    radius for an ellipsoid defined by wgs84. 
    Inputs include arrays of latitude, 
    longitude, and radius. Function returns 
    array of depths relative to this flat 
    surface with the dimensions of lon, lat. 
    """
    
    curvature = np.zeros((len(lon_Prt), \
        len(lat_Prt)), float) 
    
    return(curvature)

# %% 
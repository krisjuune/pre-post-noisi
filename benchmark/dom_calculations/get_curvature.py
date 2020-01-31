# %% Load variables
import numpy as np 
import numpy.ma as ma
import pickle 
from pathlib import Path
from math import pi
from coordinate_transformation.functions.get_domain \
    import find_nearest, truncate_domain

with open('coordinate_transformation/variables/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)

with open('coordinate_transformation/variables/lat_Prt', 'rb') as f:
    lat_Prt = pickle.load(f)
lat_Prt = np.ma.getdata(lat_Prt)

with open('coordinate_transformation/variables/bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)


# Make the domain below for calculating curvature 
# slightly bigger to have enough data (along lon)
# That way there is enough data here to cover all
lat_max = 39.5
lat_min = 35.5
lon_max = -13.8
lon_min = -19.2
bounds = [lat_max, lat_min, lon_max, lon_min]

(lat_dom, lon_dom, bathy_dom) = truncate_domain(lat_Prt, \
    lon_Prt, bathy_Prt, bounds)

# %% Get arrays with distances to the curve
import numpy as np 
from math import pi
from benchmark.dom_calculations.functions import \
    get_curvature

surface_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6370.107295)
ocean_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6365.387295)
Moho_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6357.937295)
bottom_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6270.107295)

# %% Save curvature values for the spherical case
from netCDF4 import Dataset
import numpy as np 
import datetime as dt
from coordinate_transformation.functions.get_spherical \
    import geographic_to_geocentric, wgs84
from coordinate_transformation.functions.get_rotation \
    import get_cartesian_distance
from benchmark.dom_calculations.functions import \
    get_nc_curvature

# Transform lat, lon to be centered around the N Pole
lat_N = geographic_to_geocentric(lat_Prt)
(x_N, y_N) = get_cartesian_distance(lat_N, lon_Prt)

# Save .nc datasets
get_nc_curvature('spherical_surface', surface_sphere)
get_nc_curvature('spherical_ocean', ocean_sphere)
get_nc_curvature('spherical_Moho', Moho_sphere)
get_nc_curvature('spherical_bottom', bottom_sphere)
# %% Get curvature for ellipsoid
from coordinate_transformation.transformation_functions.get_spherical \
    import radius_cnt, wgs84

def get_curvature_wgs84(lat, lon, radius = 6371, \
    theta = 38, phi = -17.75):
    """
    Function to calculate the curvature relative to 
    a flat surface at the given radius for an 
    ellipsoid defined by wgs84. Inputs include 
    arrays of latitude, longitude, and a radius. 
    Function returns array of depths relative to 
    this flat surface with the dimensions of lon, 
    lat. Units in degrees for angles, km for 
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
            # find radius at j-th latitude
            if radius == 6371: 
                a = wgs84()[0]
                b = wgs84()[1]
            else:
                # for when looking at shallower levels
                r38 = radius_cnt(theta)
                a = wgs84()[0]*radius/r38
                b = wgs84()[1]*radius/r38

            radius_j = np.sqrt((a**2*(np.cos(lat[j])**2)) + \
                (b**2*(np.sin(lat[j])**2)))/1000
            # find angle between point i,j and 
            # centre
            l1 = radius*np.sin(lon[i] - phi)
            l2 = radius*np.sin(lat[j] - theta)
            l3 = np.sqrt(np.square(l1) + \
                np.square(l2))
            # arcsin(x), x has to be [-1,1]
            alpha = np.arcsin(c/radius)
            # calculate depth to curve from flat 
            # surface
            y = radius/np.cos(alpha) - radius_j
            x = y*np.cos(alpha)
            curvature [i,j] = x 
    
    return(curvature)

surface_ellipsoid = get_curvature(lat_Prt, \
    lon_Prt, radius = 6370.107295)
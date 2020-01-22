# Laura Ermert, modified by Kristiina Joon

# Domain boundaries: 
# In [143]: lat_Prt_cnt[0]                       # Out[143]: 35.320337657545913
# In [144]: lat_Prt_cnt[-1] 
# Out[144]: 41.20709287604776
# In [145]: lat_Prt[0]                           # Out[145]: 35.502083333333331
# In [146]: lat_Prt[-1]                          # Out[146]: 41.397916666666674

#%% Transform latitudes from geographic to geocentric
import numpy as np
from math import pi, sin, cos, sqrt

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

def geographic_to_geocentric(lat):
    """
    Calculate latitude defined in the wgs84 coordinate 
    system given the geographic latitude. 
    """
    # https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    e_2 = wgs84()[2] # returns the 3rd item from function ouputs 
    lat = np.rad2deg(np.arctan((1 - e_2) * np.tan(np.deg2rad(lat))))
    return lat

#%% Express bathymetry as distance from the centre of the Earth
import numpy as np
from math import pi, sin, cos, sqrt

# Reference surface will be the Earth as defined by the WGS84 ellipsoid
def radius_cnt(lat):
    a = wgs84()[0]
    b = wgs84()[1]
    # Calculate radius for reference ellipsoid
    # in m 
    lat = pi/180*lat
    # for i in range(len(lat)): 
    r_cnt = np.sqrt((a**2*(np.cos(lat)**2)) + \
        (b**2*(np.sin(lat)**2)))
    return(r_cnt)

# The following function is basically not needed
def geocentric_to_spherical(lat_cnt, lon, elevation): 
    # Calculate radius of the reference ellipsoid
    # as a function of latitude
    m = len(lon)
    r_cnt = radius_cnt(lat_cnt)
    r_cnt = np.array([r_cnt,]*m).conj().transpose()
    # Calculate the radius for each bathymetry data point, in km 
    r_elevation = (r_cnt + elevation)/1000 
    # Calculate the colatitude to define the spherical coordinate system
    colat = 90 - lat_cnt
    return(r_elevation, colat, lon)
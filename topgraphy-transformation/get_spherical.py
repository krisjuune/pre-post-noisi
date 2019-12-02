# Laura Ermert, modified by Kristiina Joon

# Domain boundaries: 
# In [143]: lat_Prt_cnt[0]                                                                                                                                               
# Out[143]: 35.320337657545913
# In [144]: lat_Prt_cnt[-1]                                                                                                                                              
# Out[144]: 41.20709287604776
# In [145]: lat_Prt[0]                                                                                                                                                   
# Out[145]: 35.502083333333331
# In [146]: lat_Prt[-1]                                                                                                                                                  
# Out[146]: 41.397916666666674

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

def geograph_to_geocent(theta):
    # https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    e_2 = wgs84()[2] # returns the 3rd item from function ouputs 
    theta = np.rad2deg(np.arctan((1 - e_2) * np.tan(np.deg2rad(theta))))
    return theta

lat_Prt_cnt = geograph_to_geocent(lat_Prt) 

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
    r_cnt = np.sqrt((a**2*(np.cos(lat)**2)) + (b**2*(np.sin(lat)**2)))
    return(r_cnt)

len_lon = len(lon_Prt)
# Calculate radius of reference surface at geocentric latitudes
r_cnt = radius_cnt(lat_Prt_cnt) 
r_cnt = np.array([r_cnt,]*len_lon).conj().transpose()
# Calculate the radius for each bathymetry data point, in km 
r_cnt_bathy = (r_cnt + bathy_Prt)/1000 
# Calculate the colatitude to define a spherical coordinate system
colat_Prt_cnt = 90 - lat_Prt_cnt

#%% Save arrays
import numpy as np

np.save('r_cnt_bathy', r_cnt_bathy)
np.save('colat_Prt_cnt', colat_Prt_cnt)
np.save('lon_Prt', lon_Prt)

# Variables in spherical coordinates: r_cnt_bathy, colat_Prt_cnt, lon_Prt

#%%Plot spherical

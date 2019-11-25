import numpy as np
from math import pi, sin, cos, sqrt
from obspy.geodetics import gps2dist_azimuth
from pathlib import Path
from netCDF4 import Dataset

# Laura Ermert, modified by Kristiina Joon

# In [143]: lat_Prt_cnt[0]                                                                                                                                               
# Out[143]: 35.320337657545913
# In [144]: lat_Prt_cnt[-1]                                                                                                                                              
# Out[144]: 41.20709287604776
# In [145]: lat_Prt[0]                                                                                                                                                   
# Out[145]: 35.502083333333331
# In [146]: lat_Prt[-1]                                                                                                                                                  
# Out[146]: 41.397916666666674

#%%  Import global elevation files using Dataset
from pathlib import Path
from netCDF4 import Dataset

data_folder = Path('/Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/top_data/GEBCO_2019/')
file2open = data_folder / 'GEBCO_2019.nc' #file with location
nc_GEBCO = Dataset(file2open, 'r')
raw_lat = nc_GEBCO.variables['lat'][:] # import lat, in degrees N
raw_lon = nc_GEBCO.variables['lon'] [:]# import lat, in degrees E 
# import lat, in m as height above reference ellipsoid
raw_elevation = nc_GEBCO.variables['elevation'] [:]

# Define domain
#Indices for 35.5-41.4N & -22 - -14.5E obtained in MatLAB
lat_Prt = raw_lat[30120:31536] 
lon_Prt = raw_lon[37920:39720] 
bathy_Prt = raw_elevation[30120:31536, 37920:39720]
# len_lat = len(lat_Prt)
len_lon = len(lon_Prt)


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

# Calculate radius of reference surface at geocentric latitudes
r_cnt = radius_cnt(lat_Prt_cnt) 
r_cnt = np.array([r_cnt,]*len_lon).conj().transpose()
# Calculate the radius for each bathymetry data point, in m 
r_cnt_bathy = r_cnt + bathy_Prt 
# Calculate the colatitude to define a spherical coordinate system
colat_Prt_cnt = 90 - lat_Prt_cnt

#%% Save arrays
import numpy as np
import pickle
# Variables in spherical coordinates: r_cnt_bathy, colat_Prt_cnt, lon_Prt
# Load with np.load('')

r_cnt_bathy.dump('r_cnt_bathy')
colat_Prt_cnt.dump('colat_Prt_cnt')
lon_Prt.dump('lon_Prt')

#%%Plot spherical



#%% Spherical to Cartesian

def sph_to_cartesian(r,colat,lon):
    """
    Function to transform from spherical polar coordinates to 
    Cartesian coordinates. Function takes three arguments: 
    radius from origin (0,0,0), colatitude (inclination), 
    longitude (azimuth). Angles in degrees.
    """ 
    # Degrees to rad
    colat = pi/180*colat 
    lon = pi/180*lon

    # Reshape colat, lon to match shape of r
    (n,m) = r.shape

    if len(colat) == n:
        colat = np.array([colat]*m).conj().transpose()
    elif len(colat) == m: 
        colat = np.array([colat]*n)
    else:
        print('Input arguments must have at least one identical array dimension')

    if len(lon) == n:
        lon = np.array([lon]*m).conj().transpose()
    elif len(lon) == m: 
        lon = np.array([lon]*n)
    else:
        print('Input arguments must have at least one identical array dimension')

    # Transform
    x = r*np.sin(colat)*np.cos(lon)
    y = r*np.sin(colat)* np.sin(lon)
    z = r*np.cos(colat)

    return(x,y,z)


(x_Prt,y_Prt,z_Prt) = sph_to_cartesian(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

#%% Save arrays
# Variables in Cartesian coordinates: x_Prt, y_Prt, z_Prt

x_Prt.dump('x_Prt')
y_Prt.dump('y_Prt')
z_Prt.dump('z_Prt')

#%%Plot Cartesian




#%%Spherical to cylindrical
import numpy as np
from math import pi

def sph_to_cylindrical(r,colat,lon):
    """
    Function to transform from spherical polar coordinates to 
    cylindrical coordinates. Function takes three arguments: 
    radius from origin (0,0,0), colatitude (inclination), 
    longitude (azimuth). Angles in degrees.
    """ 
    # Degrees to rad
    colat = pi/180*colat 

    # Reshape colat to match shape of r
    (n,m) = r.shape

    if len(colat) == n:
        colat = np.array([colat]*m).conj().transpose()
    elif len(colat) == m: 
        colat = np.array([colat]*n)
    else:
        print('Radius and colatitude must have at least one identical array dimension')

    # Transform
    s = r*np.sin(colat)
    phi = lon
    z = r*np.cos(colat)
    
    return(s,phi,z)

(s_Prt,phi_Prt,z_Prt) = sph_to_cylindrical(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

#%% Save arrays
import numpy as np

# Variables in Cartesian coordinates: x_Prt, y_Prt, z_Prt

s_Prt.dump('s_Prt')
phi_Prt.dump('phi_Prt')
z_Prt.dump('z_Prt')

#%% Plot cylindrical 





#%% Calculate length of one degree
# # Calculate the length of one degree os each lat and lon as a function of lat

# def len_deg_lon(lat):
#     e_2 = wgs84()[2]
#     a = wgs84() [0]
#     # This is the length of one degree of longitude 
#     # approx. after WGS84, at latitude lat
#     # in m
#     lat = pi/180*lat
#     dlon = (pi*a*cos(lat))/180*sqrt((1-e_2*sin(lat)**2)) # Error only length 1 arrays can be converted to Python scalars
#     return round(dlon,5)
# def len_deg_lat(lat):
#     # This is the length of one degree of latitude 
#     # approx. after WGS84, between lat-0.5deg and lat+0.5 deg
#     # in m
#     lat = pi/180*lat
#     dlat = 111132.954 - 559.822 * cos(2*lat) + 1.175*cos(4*lat)
#     return round(dlat,5)


# dlon_Prt = len_deg_lon(lat_Prt_cnt)
# dlat_Prt = len_deg_lat(lat_Prt_cnt)

# # Why need lengths of degrees? 


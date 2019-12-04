import numpy as np
from math import pi, sin, cos, sqrt
from obspy.geodetics import gps2dist_azimuth
from pathlib import Path
import netCDF4 as nc4

# Laura Ermert, modified by Kristiina Joon
# Test for signing

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
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma

data_folder = Path('/Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/Code/coordinate_transformation/prep_data/GEBCO_2019/')
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

# Define domain
# Find indices for 35.5-41.4N & -22 - -14.5E 
lat_max = 41.4 
lat_min = 35.5
lon_max = -14.5
lon_min = -22
bounds = [lat_max, lat_min, lon_max, lon_min]

def find_nearest(array, value):
    array = np.asarray(array)
    indx = (np.abs(array - value)).argmin()
    return indx

def truncate_domain(lat, lon, value, bounds):
    """ 
    Input lat, lon, elevatin value of interest (topography, depth to Moho 
    etc.), and domain boundaries as an array of 4 elements, in the order: 
    max latitude, min latitude, max longitude, min longitude. 
    """
    #Find indices
    indx_lat_max = find_nearest(lat, bounds[0])
    indx_lat_min = find_nearest(lat, bounds[1])
    indx_lon_max = find_nearest(lon, bounds[2])
    indx_lon_min = find_nearest(lon, bounds[3])
    #Truncate domain
    lat_domain = lat[indx_lat_min:indx_lat_max]
    lon_domain = lon[indx_lon_min:indx_lon_max]
    n = len(value)
    if n == len(lat):
        value_domain = value[indx_lat_min:indx_lat_max, indx_lon_min:indx_lon_max]
    elif n == len(lon): 
        value_domain = value[indx_lon_min:indx_lon_max, indx_lat_min:indx_lat_max]
    else: 
        print('Array must have same dimensions as lat and lon')
    return(lat_domain, lon_domain, value_domain)

# Truncate domain
(lat_Prt, lon_Prt, bathy_Prt) = truncate_domain(raw_lat, raw_lon, raw_elevation, bounds)

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
    l = r*np.cos(colat)
    
    return(s,phi,l)

(s_Prt,phi_Prt,l_Prt) = sph_to_cylindrical(r_cnt_bathy, colat_Prt_cnt, lon_Prt)

#%% Save arrays
import numpy as np

# Variables in cylidrical coordinates: s_Prt, phi_Prt, z_Prt

s_Prt.dump('s_Prt')
phi_Prt.dump('phi_Prt')
z_Prt.dump('z_Prt')

#%% Plot cylindrical 





#%% Save .nc file
from pathlib import Path
from netCDF4 import Dataset
import numpy as np 
import datetime as dt

# Create .nc file
f = Dataset('topography_coord.nc','w', format='NETCDF4')
f.description = 'Togography data in sph, Cart, and cyl coordinates'

# Create dimensions
f.createDimension('colat', len(colat_Prt_cnt))
f.createDimension('lon', len(lon_Prt))

# Create variables, 'f4' for single precision floats, i.e. 32bit
radial_distance = f.createVariable('radial_distance', 'f4', ('colat', 'lon'))
colatitude = f.createVariable('colatitude', 'f4', 'colat')
longitude = f.createVariable('longitude', 'f4', 'lon')
x_value = f.createVariable('x_value', 'f4', ('colat', 'lon'))
y_value = f.createVariable('y_value', 'f4', ('colat', 'lon'))
z_value = f.createVariable('z_value', 'f4', ('colat', 'lon'))
cyl_radial_distance = f.createVariable('cyl_radial_distance', 'f4', ('colat', 'lon'))
azimuth = f.createVariable('azimuth', 'f4', 'lon')
height = f.createVariable('height', 'f4', ('colat', 'lon'))

radial_distance [:] = r_cnt_bathy
colatitude [:] = colat_Prt_cnt
longitude [:] = lon_Prt
x_value [:] = x_Prt
y_value [:] = y_Prt
z_value [:] = z_Prt
cyl_radial_distance [:] = s_Prt
azimuth [:] = phi_Prt
height [:] = l_Prt

# Add attributes to the file
today = dt.datetime.now()
f.history = "Created " + today.strftime("%d/%m/%y")
#Add local attributes to variable instances
radial_distance.units = 'km'
colatitude.units = 'degrees north'
longitude.units = 'degrees east'
x_value.units = 'km'
y_value.units = 'km'
z_value.units = 'km'
cyl_radial_distance.units = 'km'
azimuth.units = 'degrees east'
height.units = 'km'

f.close()

#%% Rotate to be centred about the N Pole
import numpy as np 
from math import pi, sin, cos
import netCDF4

#Find source lat and lon (centre of domain), calculate length of degrees
#to do it accurately but at the moment just going to take the average of 
#lat and lon as the centre

av_lat = np.mean(lat_Prt) # but need this in geographic, not geocentric
av_lon = np.mean(lon_Prt)

# Rotate
def rotation_matrix(theta,phi): 
    # Preallocate output array
    Q = np.zeros((3,3)) 

    # Fill in rotation matrix
    Q[0,0] = cos(theta)*cos(phi)
    Q[0,1] = -sin(phi)
    Q[0,2] = sin(theta)*cos(phi)
    Q[1,1] = cos(theta)*sin(phi)
    Q[1,1] = cos(phi)
    Q[1,2] = sin(theta)*sin(phi)
    Q[2,0] = -sin(theta)
    Q[2,1] = 0
    Q[2,2] = cos(theta)

    return(Q)

def rotate_N_Pole(src_lat, src_lon, x, y, z): 
    """
    Input source grographic latitude & longitude, and data file in 
    Cartesian coordinates. Rotates to N pole using rotation_matrix, 
    assuming depth of source is effectively 0. Returns the rotated 
    data in Cartesian coordinates.
    """
    # To radians
    src_lat = pi/180*src_lat
    src_lon = pi/180*src_lon

    # Get geocentric latitude, in rad
    e_2 = wgs84()[2]
    src_colat = np.arctan((1 - e_2) * np.tan(src_lat))

    # Preallocate output arrays
    (m,n) = x.shape
    x_rot = np.zeros((m,n), float)
    y_rot = np.zeros((m,n), float)
    z_rot = np.zeros((m,n), float)

    # Compute rotated x, y, z, using matrix multiplication matmul
    Q = rotation_matrix(src_colat, src_lon).transpose()
    for i in range(m): 
        a = np.concatenate((x_rot[0,], y_rot[0,], z_rot[0,]), axis = 0)
        a = a.reshape(3,n)
        b = np.concatenate((x[0,], y[0,], z[0,]), axis = 0)
        b = b.reshape(3,n)
        a = np.matmul(Q, b)
        x_rot[i,] = a[0,]
        y_rot[i,] = a[1,]
        z_rot[i,] = a[2,]
    
    return(x_rot, y_rot, z_rot)

# (x_rot, y_rot, z_rot) = rotate_N_pole(av_lat, av_lon, x_Prt, y_Prt, z_Prt)

(x_Prt_rot, y_Prt_rot, z_Prt_rot) = rotate_N_Pole(38.45, -18.2, x_Prt, y_Prt, z_Prt)

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


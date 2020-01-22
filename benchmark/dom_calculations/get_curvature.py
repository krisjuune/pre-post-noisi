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

# Make the domain be exactly 200 km in all 3 runs
# That way there is enough data here to cover all
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

surface_sphere = get_curvature(lat_Prt, \
    lon_Prt, radius = 6370.107295)
ocean_sphere = get_curvature(lat_Prt, \
    lon_Prt, radius = 6365.387295)
Moho_sphere = get_curvature(lat_Prt, \
    lon_Prt, radius = 6357.937295)
bottom_sphere = get_curvature(lat_Prt, \
    lon_Prt, radius = 6320.107295)

# %% Save curvature values for the spherical case
from netCDF4 import Dataset
import numpy as np 
import datetime as dt
from coordinate_transformation.functions.get_spherical \
    import geographic_to_geocentric, wgs84
# from coordinate_transformation.functions.get_rotation \
#     import get_cartesian_distance

def get_nc_curvature(filename, curvature_variable):
    """
    Writes a netCDF4 file with x_distance, y_distance, 
    and curvature. filename should be a string (with
    no file extension) and curvature_variable an array 
    containing the calculated curvature values. 
    """
    
    # Create .nc file
    f = Dataset(filename + '.nc','w', format = 'NETCDF4')
    f.description = 'Curvature calculated relative to tangent' + \
        ' surface at the centre of domain assuming spherical Earth'

    # Create dimensions
    f.createDimension('x', len(x_N))
    f.createDimension('y', len(y_N))

    # Create variables, 'f4' for single precision floats, i.e. 32bit
    curvature = f.createVariable('curvature', 'f4', ('x', 'y'))
    curvature [:] = curvature_variable
    x = f.createVariable('x_distance', 'f4', 'x')
    x [:] = x_N
    y = f.createVariable('y_distance', 'f4', 'y')
    y [:] = y_N

    # Add attributes to the file
    today = dt.datetime.now()
    f.history = "Created " + today.strftime("%d/%m/%y")
    #Add local attributes to variable instances
    curvature.units = 'km'
    x.units = 'km'
    y.units = 'km'

    f.close()

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
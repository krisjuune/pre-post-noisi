import numpy as np
from math import pi, sin, cos, sqrt, atan
# from coordinate_transformation.transformation_functions.get_spherical import wgs84
from coordinate_transformation.functions.get_spherical import wgs84

# %% Get distances on flat surface from source
from coordinate_transformation.functions.get_spherical \
    import radius_cnt

def get_cartesian_distance(lat, lon, \
    src_lat = 38, src_lon = -17.75):
    """
    Calculate distance of each point of lat and lon
    from the source location on a flat surface, 
    tangential to the source. Returns x (lon), y 
    (lat) in km for AxiSEMCartesian. 
    """
    # preallocate output arrays
    x = np.zeros(len(lon), float)
    y = np.zeros(len(lat), float)

    # convert to radians
    lon = pi/180*lon
    lat = pi/180*lat
    src_lon = pi/180*src_lon
    src_lat = pi/180*src_lat

    # find radius at source
    radius = radius_cnt(np.rad2deg(src_lat))/1000
    
    # find distances
    x = radius*np.sin(lon - src_lon)
    y = radius*np.sin(lat - src_lat)

    return(x,y)

#%% Define rotation matrix, Kuangdai Leng

def rotation_matrix(colat,phi): 
    """
    Colat - colatitude of source, phi - longitude of 
    source (both in radians). Function returns a 3-by-3
    rotation matrix. 
    """
    # Preallocate output array
    Q = np.zeros((3,3)) 
    # Fill in rotation matrix
    Q[0,0] = cos(colat)*cos(phi)
    Q[0,1] = -sin(phi)
    Q[0,2] = sin(colat)*cos(phi)
    Q[1,1] = cos(colat)*sin(phi)
    Q[1,1] = cos(phi)
    Q[1,2] = sin(colat)*sin(phi)
    Q[2,0] = -sin(colat)
    Q[2,1] = 0
    Q[2,2] = cos(colat)

    return(Q)

#%% Rotate source, defined by georgaphic coordinates, to N pole 

def rotate_N_pole(src_lat, src_lon, x, y, z): 
    """
    Input source grographic latitude & longitude, and data file in 
    Cartesian coordinates. Rotates to N pole using rotation_matrix. Returns the rotated 
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

#%% Rotate back to spherical coordinates for plotting

def cartesian_to_sph(x,y,z):
    """
    Input rotated data in Cartesian coordinates, transform to 
    spherical to input into AxiSEM. 
    Return r (in km), colatitude and longitude (in degrees).
    """
    # To spherical 
    r = sqrt(x^2 + y^2 + z^2)
    colat = arctan(sqrt(x^2 + y^2)/z)
    phi = arctan(y/x)
    # To degrees
    colat = 180/pi*colat
    phi = 180/pi*phi
    return(r,colat,phi)

#%% Get elevation relative to reference (average) elevation
# The input file for AxiSEM needs the elevation data for 
# relabelling relative to a reference depth, down is positive, 
# up is negative. 

def rel_depth(r):
    """
    Given the data as distance from the centre of the Earth, 
    return depths relative to the 'average depth'. 
    """
    # Do I need to account for ellipticity again? 
    z = r # z is some function of r
    return(z)
import numpy as np
from math import pi, sin, cos, sqrt
# from coordinate_transformation.transformation_functions.get_spherical import wgs84
from transformation_functions.get_spherical import wgs84

#%% Define rotation matrix, Kuangdai Leng

def rotation_matrix(theta,phi): 
    """
    Theta - colatitude of source, phi - longitude of 
    source (both in radians). Function returns a 3-by-3
    rotation matrix. 
    """
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

#%% Rotate source, defined by georgaphic coordinates, to N pole 

def rotate_N_pole(src_lat, src_lon, x, y, z): 
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

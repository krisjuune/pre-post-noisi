import numpy as np 
from math import pi, sin, cos
import netCDF4

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


#%% Calculate colatitude of source
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

def lat_to_colat(lat, depth, geocentric): 
    """
    Due to ellipticity, the colatitude cannot just be calculated by 
    the simple 90-lat method. Function returns colatitude, given 
    either the geographic latitude or the geocetric latitude, and 
    depth of source. Set geocentric equal to 'True' if input latitudes 
    are geocentric, and 'False' if geographic. 
    """
    if geocentric != 'True' or 'False':
        print('geocentric has to be a Boolean value')
    elif geocentric == 'True':
        # Ellipticity already accounted for
        if depth == 0:
            colat = 90 - lat
        else: 
            print('Function only works for 0 depths')
    # else:
    #     e_2 = wgs84()[2] # returns the 3rd item from function ouputs 
    #     if depth == 0: 
    #         colat = np.rad2deg(np.arctan((1 - e_2) * np.tan(np.deg2rad(lat))))
    #         # This is for source depth = 0, how to incorporate depths?
    #     else:
    #         print('Function only works for 0 depths')

    return(colat)


#%% Rotate to N Pole
import numpy as np
from math import pi, sin, cos, sqrt

# From Kuangdai's code: 
# RDCol3 Geodesy::rotateGlob2Src(const RDCol3 &rtpG, double srclat, double srclon, double srcdep) {
#     const RDCol3 &xyzG = toCartesian(rtpG);
#     const RDCol3 &xyzS = rotationMatrix(lat2Theta_d(srclat, srcdep), lon2Phi(srclon)).transpose() * xyzG;
#     bool defined = true;
#     RDCol3 rtpS = toSpherical(xyzS, defined);
#     if (!defined) {
#         rtpS(2) = rtpG(2);
#     }
#     return rtpS;
# }

### SCHEISSE CODE
# def rotate_N_pole(src_lat, src_lon, src_depth, x, y, z):
#     """
#     Input source latitude & longitude, and data file in Cartesian 
#     coordinates. Rotates to N pole using the rotation matrix. Returns
#     the rotated data in Cartesian coordinates.
#     """

#     # Calculate colatitude of source given depth and geographic lat
#     src_colat = lat_to_colat(src_lat, src_depth, 'False')
#     # Longitude in radians
#     src_lon = pi/180*src_lon

#     # How would this multipliciation work as rotation matrix is
#     # 3-by-3 but each of x, y, and z is len(lat)-by-len(lon)
#     x_rot = rotation_matrix(src_colat, src_lon).transpose() * x
#     y_rot = rotation_matrix(src_colat, src_lon).transpose() * y
#     z_rot = rotation_matrix(src_colat, src_lon).transpose() * z

#     return(x_rot, y_rot, z_rot)

# # Get source centered in Cartesian

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

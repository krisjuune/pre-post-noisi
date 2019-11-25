import numpy as np 
from math import pi, sin, cos
import netCDF4
import netcdf

#%% Define rotation matrix

# From Kuangdai's code: 
# RDMat33 Geodesy::rotationMatrix(double theta, double phi) {
#     RDMat33 Q;
#     Q(0, 0) = cos(theta) * cos(phi);
#     Q(0, 1) = -sin(phi);
#     Q(0, 2) = sin(theta) * cos(phi);
#     Q(1, 0) = cos(theta) * sin(phi);
#     Q(1, 1) = cos(phi);
#     Q(1, 2) = sin(theta) * sin(phi);
#     Q(2, 0) = -sin(theta);
#     Q(2, 1) = 0.;
#     Q(2, 2) = cos(theta);
#     return Q;
# }

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

#%% Rotate to N Pole

# From Kunagdai's code: 
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

def rotate_NPole(src_lat, src_lon, x, y, z):
    x_rot = rotation_matrix(src_lat, src_lon).transpose() * x
    y_rot = rotation_matrix(src_lat, src_lon).transpose() * y
    z_rot = rotation_matrix(src_lat, src_lon).transpose() * z
    return(x_rot, y_rot, z_rot)

# Need theta and phi of the source in the rotation matrix, 
# then multiply the rot matrix result transpose wiht Cart
# coordinates of the model 

# Get source centered in Cartesian
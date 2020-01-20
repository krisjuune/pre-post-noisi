import numpy as np
from math import sin, cos, pi

#%% Transform into cylindrical 

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

#%% Plot cylindrical

import numpy as np
from math import sin, cos, pi

# # Load saved input files
# r_cnt_bathy_ls = np.load('r_cnt_bathy.npy')
# colat_Prt_cnt_ls = np.load('colat_Prt_cnt.npy')
# lon_Prt_ls = np.load('lon_Prt.npy')

# # Convert from np.str to np.float
# r_cnt_bathy = r_cnt_bathy_ls.astype(np.float)
# colat_Prt_cnt = colat_Prt_cnt_ls.astype(np.float)
# lon_Prt = lon_Prt_ls.astype(np.float)

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

# %% Plot Cartesian

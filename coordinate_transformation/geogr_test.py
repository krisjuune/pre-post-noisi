import numpy as np
from math import pi, sin, cos, sqrt
from obspy.geodetics import gps2dist_azimuth

# Laura Ermert, modified by Kristiina Joon

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
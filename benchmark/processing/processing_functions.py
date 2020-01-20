import numpy as np 
from pathlib import Path
from math import pi
from coordinate_transformation.transformation_functions.get_spherical \
    import wgs84

def station_data(path, station):
    """
    Function that retrieves the seismic data from station
    'station' given the relative path 'path', both inputs 
    are strings. This works for II type (not IU) stations.
    Function returns the data array. 
    """
    path = Path(path)
    file = 'II.' + station + '.RTZ.ascii'
    # file handle
    file = path/file 
    # Open file and retrieve data
    raw_data = open(file, 'r')
    raw_data = raw_data.read()
    raw_data = np.fromstring(raw_data, dtype = float, sep=' ')
    # nr of columns is always 4 since time, rr, tt, zz
    m = int(len(raw_data)/4) 
    # preallocate output array
    data = np.zeros(((m),4), float)
    # retrieve data (which has been sorted row-wise)
    # and sort it column-wise, returning every 4th element
    data[:,0] = raw_data[0::4] 
    data[:,1] = raw_data[1::4]
    data[:,2] = raw_data[2::4]
    data[:,3] = raw_data[3::4]
    return(data)

# Calculate the length of one degree of lat and lon as a function of lat
def len_deg_lon(lat):
    """
    Calculates length of one degree of longitude
    at latitudes lat. Input lat must be an array
    of integers. 
    """
    e_2 = wgs84()[2]
    a = wgs84() [0]
    # This is the length of one degree of longitude 
    # approx. after WGS84, at latitude lat
    # in m
    lat = pi/180*lat
    dlon = (pi*a*np.cos(lat))/180*np.sqrt((1-e_2*np.sin(lat)**2))
    return np.round(dlon,5)

def len_deg_lat(lat):
    """
    Calculates length of one degree of latitude
    at latitudes lat. Input lat must be an array
    of integers. 
    """
    # This is the length of one degree of latitude 
    # approx. after WGS84, between lat-0.5deg and lat+0.5 deg
    # in m
    lat = pi/180*lat
    dlat = 111132.954 - 559.822 * np.cos(2*lat) + 1.175*np.cos(4*lat)
    return np.round(dlat,5)

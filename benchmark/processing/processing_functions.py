import numpy as np 
from pathlib import Path

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


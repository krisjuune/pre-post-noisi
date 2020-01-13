import numpy as np 
import obspy as obs 
from pathlib import Path
import matplotlib.pyplot as plt

# %% Get the data
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

path1 = 'validation/processing/raw_data/sim1/'
path2 = 'validation/processing/raw_data/sim2/'
path3 = 'validation/processing/raw_data/sim3/'
station = 'LAT5'
data1 = station_data(path1, station)
data2 = station_data(path2, station)
data3 = station_data(path3, station)

station_lon = 'LON5'
data1_lon = station_data(path1, station_lon)
data2_lon = station_data(path2, station_lon)
data3_lon = station_data(path3, station_lon)
# %% Plot seismograms
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

plt.figure(1)
fig = plt.subplot(211)
plt.plot(data1[0:400,0], data1[0:400,3], color = 'orange', \
    linestyle = '--', linewidth = '1')
plt.plot(data2[0:400,0], data2[0:400,3], color = 'orange', \
    linewidth = '1')
plt.plot(data3[0:400,0], data3[0:400,3], color = 'darkblue', \
    linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('Lat')
fig.axes.get_xaxis().set_visible(False)

plt.subplot(212)
plt.plot(data1_lon[0:400,0], data1_lon[0:400,3], color = 'orange', \
    linestyle = '--', linewidth = '1')
plt.plot(data2_lon[0:400,0], data2_lon[0:400,3], color = 'orange', \
    linewidth = '1')
plt.plot(data3_lon[0:400,0], data3_lon[0:400,3], color = 'darkblue', \
    linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('Lon')

plt.savefig('validation_test.png', dpi = 600)

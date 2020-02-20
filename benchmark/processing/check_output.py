# %%
import numpy as np 
# import obspy as obs 
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

path1 = 'benchmark/processing/raw_data/cartesian/'
path2 = 'benchmark/processing/raw_data/spherical/'
path3 = 'benchmark/processing/raw_data/geographic/'
station = 'ST3'
data1 = station_data(path1, station)
data2 = station_data(path2, station)
data3 = station_data(path3, station)

station_lon = 'ST1'
data1_lon = station_data(path1, station_lon)
data2_lon = station_data(path2, station_lon)
data3_lon = station_data(path3, station_lon)

station_w = 'ST0'
data1_w = station_data(path1, station_w)
data2_w = station_data(path2, station_w)
data3_w = station_data(path3, station_w)

# %% Plot seismograms
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

# TODO look at the code for plotting curvature
# it has better handling of axes stuff so maybe need to 
# set some variable to axes of each subplot separately

plt.figure(1)
# set desired range of data (time plotted)
m = np.arange(0,4000)
# set desired compenent, 1 radial, 2 transverse, 3 vertical
n = 3
lwidth = '0.5'

fig = plt.subplot(311)
plt.plot(data1[:,0], data1[:,n]*(-1), color = 'orange', \
    linestyle = '--', linewidth = lwidth)
plt.plot(data2[:,0], data2[:,n], color = 'orange', \
    linewidth = lwidth)
plt.plot(data3[:,0], data3[:,n], color = 'darkblue', \
    linewidth = lwidth, alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('170 km', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)

fig1 = plt.subplot(312)
plt.plot(data1_lon[:,0], data1_lon[:,n]*(-1), color = 'orange', \
    linestyle = '--', linewidth = lwidth)
plt.plot(data2_lon[:,0], data2_lon[:,n], color = 'orange', \
    linewidth = lwidth)
plt.plot(data3_lon[:,0], data3_lon[:,n], color = 'darkblue', \
    linewidth = lwidth, alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('60 km', fontsize = 10)
fig1.axes.get_xaxis().set_visible(False)

plt.subplot(313)
plt.plot(data1_w[:,0], data1_w[:,n]*(-1), color = 'orange', \
    linestyle = '--', linewidth = lwidth)
plt.plot(data2_w[:,0], data2_w[:,n], color = 'orange', \
    linewidth = lwidth)
plt.plot(data3_w[:,0], data3_w[:,n], color = 'darkblue', \
    linewidth = lwidth, alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('at source', fontsize = 10)
plt

plt.subplots_adjust(hspace=0.4)
plt.show()
plt.savefig('which_curvature.png', dpi = 600)

# %%
fig, axs = plt.subplots(3, sharex=True, sharey=True, gridspec_kw={'hspace': 0})

fig.suptitle('Sharing both axes')
axs[0].plot(x, y ** 2)
axs[1].plot(x, 0.3 * y, 'o')
axs[2].plot(x, y, '+')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
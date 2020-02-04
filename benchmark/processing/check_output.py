#%%
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
# path3 = 'Cartesian_3layers/output/stations/'
station = 'ST3'
data1 = station_data(path1, station)
data2 = station_data(path2, station)
# data3 = station_data(path3, 'ST1')

station_lon = 'ST1'
data1_lon = station_data(path1, station_lon)
data2_lon = station_data(path2, station_lon)
# data3_lon = station_data(path3, station_lon)

station_w = 'ST0'
data1_w = station_data(path1, station_w)
data2_w = station_data(path2, station_w)
# data3_lon = station_data(path3, station_lon)
# %% Plot seismograms
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

plt.figure(1)
# set desired range of data (time plotted)
m = np.arange(0,4000)
# set desired compenent, 1 radial, 2 transverse, 3 vertical
n = 3
fig = plt.subplot(311)
plt.plot(data1[:,0], data1[:,3]*(-1), color = 'orange', \
    linestyle = '--', linewidth = '1')
plt.plot(data2[:,0], data2[:,n], color = 'orange', \
    linewidth = '1')
# plt.plot(data3[m,0], data3[m,n], color = 'darkblue', \
#     linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('170 km', fontsize = 10)
# fig.axes.get_xaxis().set_visible(False)

plt.subplot(312)
plt.plot(data1_lon[:,0], data1_lon[:,3]*(-1), color = 'orange', \
    linestyle = '--', linewidth = '1')
plt.plot(data2_lon[:,0], data2_lon[:,n], color = 'orange', \
    linewidth = '1')
# plt.plot(data3_lon[m,0], data3_lon[m,n], color = 'darkblue', \
#     linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('60 km', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_xaxis().set_visible(False)

plt.subplot(313)
plt.plot(data1_w[:,0], data1_w[:,3]*(-1), color = 'orange', \
    linestyle = '--', linewidth = '1')
plt.plot(data2_w[:,0], data2_w[:,n], color = 'orange', \
    linewidth = '1')
# plt.plot(data3_lon[m,0], data3_lon[m,n], color = 'darkblue', \
#     linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('at source', fontsize = 10)
plt

plt.subplots_adjust(hspace=0.4)
plt.show()
# plt.savefig('full_domain.png', dpi = 600)

# %%
fig, axs = plt.subplots(3, sharex=True, sharey=True, gridspec_kw={'hspace': 0})

fig.suptitle('Sharing both axes')
axs[0].plot(x, y ** 2)
axs[1].plot(x, 0.3 * y, 'o')
axs[2].plot(x, y, '+')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
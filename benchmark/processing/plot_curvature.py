import numpy as np 
# import obspy as obs 
from pathlib import Path
import matplotlib.pyplot as plt
# plt.style.use('ggplot')

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

path1 = 'benchmark/processing/raw_data/curvature/cartesian/stations'
path2 = 'benchmark/processing/raw_data/curvature/spherical/stations'
path3 = 'benchmark/processing/raw_data/curvature/geographic/stations'
station = 'ST4' 
data1 = station_data(path1, station)
data2 = station_data(path2, station)
data3 = station_data(path3, station)

station_lon = 'ST2'
data1_lon = station_data(path1, station_lon)
data2_lon = station_data(path2, station_lon)
data3_lon = station_data(path3, station_lon)

station_w = 'ST0'
data1_w = station_data(path1, station_w)
data2_w = station_data(path2, station_w)
data3_w = station_data(path3, station_w)

# 210, 170, 130, 85, 40 km for the station distances

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
lwidth = 0.75

fig = plt.subplot(311)
plt.plot(data1_w[:,0], data1_w[:,n]*(-1), color = 'orange', \
    linestyle = '--', linewidth = lwidth)
plt.plot(data2_w[:,0], data2_w[:,n], color = 'orange', \
    linewidth = lwidth)
plt.plot(data3_w[:,0], data3_w[:,n], color = 'darkblue', \
    linewidth = lwidth, alpha = 0.75)
# lgd = fig.axes.legend(('Cartesian', 'Spherical', 'Geographic'), \
#     prop={'size': 6}, bbox_to_anchor=(1, 1.05), frameon=False)
# lgd.get_frame().set_edgecolor('k')
plt.ylim(-199999, 199999)
plt.ylabel('Displ.')
plt.title('0 km', fontsize = 10)
fig.axes.tick_params(labelbottom=False, bottom=False, 
                     labelleft=False, left=False) 

fig1 = plt.subplot(312)
plt.plot(data1_lon[:,0], data1_lon[:,n]*(-1), color = 'orange', \
    linestyle = '--', linewidth = lwidth)
plt.plot(data2_lon[:,0], data2_lon[:,n], color = 'orange', \
    linewidth = lwidth)
plt.plot(data3_lon[:,0], data3_lon[:,n], color = 'darkblue', \
    linewidth = lwidth, alpha = 0.75)
lgd1 = fig1.axes.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 8}, bbox_to_anchor=(1, 0.96), frameon=False)
# lgd1.get_frame().set_edgecolor('k')
plt.ylim(-3900, 3900)
plt.ylabel('Displ.')
plt.title('85 km', fontsize = 10)
fig1.axes.tick_params(labelbottom=False, bottom=False, 
                      labelleft=False, left=False) 

fig2 = plt.subplot(313)
plt.plot(data1[:,0], data1[:,n]*(-1), color = 'orange', \
    linestyle = '--', linewidth = lwidth)
plt.plot(data2[:,0], data2[:,n], color = 'orange', \
    linewidth = lwidth)
plt.plot(data3[:,0], data3[:,n], color = 'darkblue', \
    linewidth = lwidth, alpha = 0.75)
# lgd2 = fig2.axes.legend(('Cartesian', 'Spherical', 'Geographic'), \
#     prop={'size': 6}, bbox_to_anchor=(1, 1.05))
# lgd2.get_frame().set_edgecolor('k')
plt.ylim(-1300,1300)
plt.ylabel('Displ.')
plt.title('170 km', fontsize = 10)
plt.xlabel('Time (s)')
fig2.axes.tick_params(labelleft=False, left=False) 

plt.subplots_adjust(hspace=0.35)
plt.savefig('which_curvature.png', dpi = 900, bbox_extra_artists=(lgd1,), bbox_inches='tight')
plt.show()

# %%

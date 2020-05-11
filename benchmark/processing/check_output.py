# %%
import numpy as np 
# import obspy as obs 
from pathlib import Path
import matplotlib.pyplot as plt
from netCDF4 import Dataset

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

def synthetics_data(path, station):
    """
    Load one station data from synthetics.nc. 
    """
    path = Path(path)
    nc = Dataset(path/'axisem3d_synthetics.nc')
    data = np.asarray(nc.variables['II.'+station+'.ENZ'][:,:])
    time = np.asarray(nc.variables['time_points'][:])
    return(data, time)

path1 = '../../../runs/benchmark/test_curvature/' + 'Cartesian_3layers' + '/output/stations/'
path2 = '../../../runs/benchmark/test_curvature/' + 'Spherical_3layers' + '/output/stations/'
path3 = '../../../runs/benchmark/test_curvature/' + 'Geographic_3layers' + '/output/stations/'

station_mid = 'ST4'
data1_mid, time1 = synthetics_data(path1, station_mid)
data2_mid, time2 = synthetics_data(path2, station_mid)
data3_mid, time3 = synthetics_data(path3, station_mid)

station_far = 'ST2'
data1_far, time1 = synthetics_data(path1, station_far)
data2_far, time2 = synthetics_data(path2, station_far)
data3_far, time3 = synthetics_data(path3, station_far)

station_src = 'ST0'
data1_src, time1 = synthetics_data(path1, station_src)
data2_src, time2 = synthetics_data(path2, station_src)
data3_src, time3 = synthetics_data(path3, station_src)

# %% Plot seismograms
import matplotlib.pyplot as plt
import numpy as np 

# TODO look at the code for plotting curvature
# it has better handling of axes stuff so maybe need to 
# set some variable to axes of each subplot separately

plt.figure(1)
# set desired range of data (time plotted)
m = -1
# set desired compenent, 0 radial, 1 transverse, 2 vertical (or ENZ)
n = 2
lwidth = '0.5'

fig = plt.subplot(311)
plt.plot(time1[:-2570], data1_mid[:-2570,n]*(-1), color = 'orange', linestyle = '--', linewidth = lwidth)
plt.plot(time2[:-2570], data2_mid[:-2570,n], color = 'orange', linewidth = lwidth)
plt.plot(time3[:-2570], data3_mid[:-2570,n], color = 'darkblue', linewidth = lwidth, alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('170 km', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)

fig1 = plt.subplot(312)
plt.plot(time1[:-2570], data1_far[:-2570,n]*(-1), color = 'orange', linestyle = '--', linewidth = lwidth)
plt.plot(time2[:-2570], data2_far[:-2570,n], color = 'orange', linewidth = lwidth)
plt.plot(time3[:-2570], data3_far[:-2570,n], color = 'darkblue', linewidth = lwidth, alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('60 km', fontsize = 10)
fig1.axes.get_xaxis().set_visible(False)

plt.subplot(313)
plt.plot(time1[:-2570], data1_src[:-2570,n]*(-1), color = 'orange', linestyle = '--', linewidth = lwidth)
plt.plot(time2[:-2570], data2_src[:-2570,n], color = 'orange', linewidth = lwidth)
plt.plot(time3[:-2570], data3_src[:-2570,n], color = 'darkblue', linewidth = lwidth, alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('at source', fontsize = 10)
plt

plt.subplots_adjust(hspace=0.4)
plt.savefig('which_curvature.png', dpi = 600)
plt.show()

# %% Plot reg vs oversampled TODO look completely different, should be the same
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.filter import *

path4 = '../../../runs/benchmark/test_curvature/' + 'Geographic_2.5mesh' + '/output/stations/'
data4_mid, time4 = synthetics_data(path4, station_mid)
data4_far, time4 = synthetics_data(path4, station_far)
data4_src, time4 = synthetics_data(path4, station_src)

# path5 = '../../../runs/benchmark/test_curvature/' + 'Geographic_1mesh' + '/output/stations/'
# data5_mid, time5 = synthetics_data(path5, station_mid)
# data5_far, time5 = synthetics_data(path5, station_far)
# data5_src, time5 = synthetics_data(path5, station_src)

plt.figure(1)
n = 2
lwidth = '0.5'

fig = plt.subplot(311)
# plt.plot(time5, data5_mid[:,n], color = 'skyblue', linewidth = lwidth)
plt.plot(time4, data4_mid[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(time3, data3_mid[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('2.5s mesh', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('170 km', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)

fig1 = plt.subplot(312)
# plt.plot(time5, data5_far[:,n], color = 'skyblue', linewidth = lwidth)
plt.plot(time4, data4_far[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(time3, data3_far[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('2.5s mesh', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('60 km', fontsize = 10)
fig1.axes.get_xaxis().set_visible(False)

plt.subplot(313)
# plt.plot(time5, data5_src[:,n], color = 'skyblue', linewidth = lwidth)
plt.plot(time4, data4_src[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(time3, data3_src[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('below source', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('210 km', fontsize = 10)
plt

plt.subplots_adjust(hspace=0.4)
plt.savefig('oversampled.png', dpi = 600)
plt.show()


# %%

import numpy as np 
import obspy as obs 
from pathlib import Path
import matplotlib.pyplot as plt

# %% Get the data
from benchmark.processing.processing_functions \
    import station_data

path1 = 'benchmark/processing/raw_data/cartesian/'
path2 = 'benchmark/processing/raw_data/spherical/'
# path3 = 'raw_data/sim3/'
station = 'LAT4'
data1 = station_data(path1, station)
data2 = station_data(path2, station)
# data3 = station_data(path3, station)

station_lon = 'LON4'
data1_lon = station_data(path1, station_lon)
data2_lon = station_data(path2, station_lon)
# data3_lon = station_data(path3, station_lon)
# %% Plot seismograms
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

plt.figure(1)
fig = plt.subplot(211)
plt.plot(data1[0:1000,0], data1[0:1000,3], color = 'orange', \
    linestyle = '-', linewidth = '1')
plt.plot(data2[0:1000,0], data2[0:1000,3], color = 'darkblue', \
    linewidth = '1', alpha = 0.65)
# plt.plot(data3[:,0], data3[:,3], color = 'darkblue', \
    # linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('Lat')
fig.axes.get_xaxis().set_visible(False)

plt.subplot(212)
plt.plot(data1_lon[0:1000,0], data1_lon[0:1000,3], color = 'orange', \
    linestyle = '-', linewidth = '1')
plt.plot(data2_lon[0:1000,0], data2_lon[0:1000,3], color = 'darkblue', \
    linewidth = '1', alpha = 0.65)
# plt.plot(data3_lon[:,0], data3_lon[:,3], color = 'darkblue', \
#     linewidth = '1', alpha = 0.65)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6})
plt.title('Lon')

# plt.savefig('validation_test.png', dpi = 600)


# %%

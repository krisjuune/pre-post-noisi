import numpy as np 
import obspy as obs 
from pathlib import Path
import matplotlib.pyplot as plt

# %% Get the data
from benchmark.processing.processing_functions \
    import station_data

path1 = 'raw_data/sim1.1/output/stations'
path2 = 'raw_data/sim2/'
path3 = 'raw_data/sim3/'
station = 'LAT6'
data1 = station_data(path1, station)
data2 = station_data(path2, station)
data3 = station_data(path3, station)

station_lon = 'LAT2'
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
# plt.plot(data2[0:400,0], data2[0:400,3], color = 'orange', \
#     linewidth = '1')
# plt.plot(data3[0:400,0], data3[0:400,3], color = 'darkblue', \
#     linewidth = '1', alpha = 0.65)
# plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
#     prop={'size': 6})
plt.title('far')
fig.axes.get_xaxis().set_visible(True)

plt.subplot(212)
plt.plot(data1_lon[0:400,0], data1_lon[0:400,3], color = 'orange', \
    linestyle = '--', linewidth = '1')
# plt.plot(data2_lon[0:400,0], data2_lon[0:400,3], color = 'orange', \
#     linewidth = '1')
# plt.plot(data3_lon[0:400,0], data3_lon[0:400,3], color = 'darkblue', \
#     linewidth = '1', alpha = 0.65)
# plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
#     prop={'size': 6})
plt.title('near')

# plt.savefig('validation_test.png', dpi = 600)


# %%

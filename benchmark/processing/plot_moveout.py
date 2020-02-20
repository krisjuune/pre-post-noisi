# %% Moveout plot
# Plot seismograms on time vs distance

import numpy as np 
import obspy # allows to filter as lowpass
from pathlib import Path
from benchmark.functions import station_data
import matplotlib.pyplot as plt 
    
station_45 = ('ST1', 'ST2', 'ST3', 'ST4', 'ST5')
station_135 = ('ST6', 'ST7', 'ST8', 'ST9', 'ST10')
station_225 = ('ST11', 'ST12', 'ST13', 'ST14', 'ST15')
station_315 = ('ST16', 'ST17', 'ST18', 'ST19', 'ST20')
path = 'benchmark/processing/raw_data/sensitivity/5sec'

for i in np.arange(len(station_45)):
    data_45 = station_data(path, station_45 [i])
    data_135 = station_data(path, station_135 [i])
    data_225 = station_data(path, station_225 [i])
    data_315 = station_data(path, station_315 [i])

    # choose component to be plotted
    n = 3
    lwidth = 0.75
    if i == 0: 
        plt.plot(data_45[:,0], data_45[:,n], linewidth = lwidth, color = 'orange')
        plt.plot(data_135[:,0], data_135[:,n], linewidth = lwidth, color = 'olive')
        plt.plot(data_225[:,0], data_225[:,n], linewidth = lwidth, color = 'skyblue', alpha = 0.7)
        plt.plot(data_315[:,0], data_315[:,n], linewidth = lwidth, color = 'darkblue', alpha = 0.7, linestyle = '--')
        dx = 1.8*max(data_45[:,n])
        axes = plt.gca()
        axes.set_ylim([-1.2*max(data_45[:,n]),(len(station_135)-0.3)*dx])
    else: 
        plt.plot(data_45[:,0], data_45[:,n] + i*dx, linewidth = lwidth, color = 'orange')
        plt.plot(data_135[:,0], data_135[:,n] + i*dx, linewidth = lwidth, color = 'olive')
        plt.plot(data_225[:,0], data_225[:,n] + i*dx, linewidth = lwidth, color = 'skyblue', alpha = 0.7)
        plt.plot(data_315[:,0], data_315[:,n] + i*dx, linewidth = lwidth, color = 'darkblue', alpha = 0.7, linestyle = '--')

plt.legend(('NE', 'NW', 'SW', 'SE'), prop={'size': 8}, loc='upper left')

# %%

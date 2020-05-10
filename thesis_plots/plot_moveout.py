# %% Moveout plot for 5s with bathymetry and 5s mesh
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
ticks = np.zeros(len(station_45))
labels = ['40', '85', '130', '170', '210']

for i in np.arange(len(station_45)):
    data_45 = station_data(path, station_45 [i])
    data_135 = station_data(path, station_135 [i])
    data_225 = station_data(path, station_225 [i])
    data_315 = station_data(path, station_315 [i])

    # choose component to be plotted
    n = 3
    lwidth = 0.75
    c1 = 'orange'
    c2 = 'orange'
    c3 = 'orange'
    c4 = 'darkblue'
    # c1 = 'black'
    # c2 = 'grey'
    # c3 = 'black'
    # c4 = 'grey'
    l1 = 'solid'
    l2 = 'dotted'
    l3 = 'dashdot'
    l4 = 'solid'

    if i == 0: 
        plt.plot(data_45[:,0], data_45[:,n], linewidth = lwidth, color = c1, linestyle = l1)
        plt.plot(data_135[:,0], data_135[:,n], linewidth = lwidth, color = c2, linestyle = l2, alpha = 1)
        plt.plot(data_225[:,0], data_225[:,n], linewidth = lwidth, color = c3, linestyle = l3, alpha = 1)
        plt.plot(data_315[:,0], data_315[:,n], linewidth = lwidth, color = c4, linestyle = l4, alpha = 1)
        dx = 2.0*max(data_135[:,n])
        axes = plt.gca()
        axes.set_ylim([-1.2*max(data_45[:,n]),(len(station_135)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
    else: 
        plt.plot(data_45[:,0], data_45[:,n] + i*dx, linewidth = lwidth, color = c1, linestyle = l1)
        plt.plot(data_135[:,0], data_135[:,n] + i*dx, linewidth = lwidth, color = c2, linestyle = l2, alpha = 0.90)
        plt.plot(data_225[:,0], data_225[:,n] + i*dx, linewidth = lwidth, color = c3, linestyle = l3, alpha = 0.85)
        plt.plot(data_315[:,0], data_315[:,n] + i*dx, linewidth = lwidth, color = c4, linestyle = l4, alpha = 0.90)
    ticks [i] = dx*i

plt.legend(('NE', 'NW', 'SW', 'SE'), prop={'size': 8}, loc='upper left')
axes.set_yticks(ticks)
axes.set_yticklabels(labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
axes.set_ylabel('Distance from sourcs (km)')

plt.savefig('moveout_quadrants.png', dpi = 600)

# %% Moveout plot for src sensitivity to bathymetry on 5s mesh
# wanna see convergence to flat data as source period increases
import numpy as np 
import obspy # allows to filter as lowpass
from pathlib import Path
from benchmark.functions import station_data
import matplotlib.pyplot as plt 

# TODO add comparison to flat which should be the easiest to read one on graph
# at the moment change ST10 to 16 and ST12 to 19 in geographic, doesnt work well, rerun or smth
station_315 = ('ST16', 'ST17', 'ST18', 'ST19', 'ST20')
# path1 = 'benchmark/processing/raw_data/sensitivity/25sec'
# BUG in data, 2.5 messes it up
# path0 = 'benchmark/processing/raw_data/sensitivity/25sec'
path1 = 'benchmark/processing/raw_data/sensitivity/5sec'
path2 = 'benchmark/processing/raw_data/sensitivity/10sec'
path3 = 'benchmark/processing/raw_data/sensitivity/20sec'
path4 = 'benchmark/processing/raw_data/curvature/geographic/stations'
station_315 = ('ST16', 'ST17', 'ST18', 'ST19', 'ST20')
ticks = np.zeros(len(station_315))
labels = ['40', '85', '130', '170', '210']

for i in np.arange(len(station_315)):
    # data25 = station_data(path0, station_315 [i])
    data_5 = station_data(path1, station_315 [i])
    data_10 = station_data(path2, station_315 [i])
    data_20 = station_data(path3, station_315 [i])
    data_flat = station_data(path4, station_315 [i])

    # choose component to be plotted
    n = 3
    lwidth = 0.75
    if i == 0: 
        # plt.plot(data25[:,0], data25[:,n], linewidth = lwidth, color = 'orange')
        plt.plot(data_flat[:,0], data_flat[:,n], linewidth = lwidth, color = 'dimgrey', alpha = 0.90)
        plt.plot(data_5[:,0], data_5[:,n], linewidth = lwidth, color = 'royalblue', alpha = 0.90)
        plt.plot(data_10[:,0], data_10[:,n], linewidth = lwidth, color = 'orange', alpha = 0.95)
        # plt.plot(data_20[:,0], data_20[:,n], linewidth = lwidth, color = 'royalblue', alpha = 0.65)
        dx = 1.8*max(data_5[:,n])
        axes = plt.gca()
        axes.set_ylim([-1.2*max(data_5[:,n]),(len(station_315)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
    else: 
        # plt.plot(data25[:,0], data25[:,n] + i*dx, linewidth = lwidth, color = 'orange')
        plt.plot(data_flat[:,0], data_flat[:,n] + i*dx, linewidth = lwidth, color = 'dimgrey', alpha = 0.90)
        plt.plot(data_5[:,0], data_5[:,n] + i*dx, linewidth = lwidth, color = 'royalblue', alpha = 0.90)
        plt.plot(data_10[:,0], data_10[:,n] + i*dx, linewidth = lwidth, color = 'orange', alpha = 0.95)
        # plt.plot(data_20[:,0], data_20[:,n] + i*dx, linewidth = lwidth, color = 'royalblue', alpha = 0.65)
    ticks [i] = dx*i

plt.legend(('flat', '5s', '10s'), prop={'size': 8}, loc='upper left')
axes.set_yticks(ticks)
axes.set_yticklabels(labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
axes.set_ylabel('Distance from sourcs (km)')

# plt.savefig('test.png', dpi = 600)

# %% Plot moveout for geographic flat oversampled
import numpy as np 
import obspy # allows to filter as lowpass
from pathlib import Path
from benchmark.functions import station_data
import matplotlib.pyplot as plt 

station_315 = ('ST16', 'ST17', 'ST18', 'ST19', 'ST20')
path = 'benchmark/processing/raw_data/curvature/geographic_oversampled/stations/'
ticks = np.zeros(len(station_315))
labels = ['40', '85', '130', '170', '210']

for i in np.arange(len(station_315)):
    data_315 = station_data(path, station_315 [i])

    # choose component to be plotted
    n = 3
    lwidth = 0.75
    # c1 = 'dimgrey'
    # c2 = 'darkblue'
    # c3 = 'skyblue'
    # c4 = 'orange'
    c1 = 'black'
    # c2 = 'grey'
    # c3 = 'black'
    # c4 = 'grey'
    l1 = 'solid'
    # l2 = 'solid'
    # l3 = 'dashed'
    # l4 = 'dashed'

    if i == 0: 
        plt.plot(data_315[:,0], data_315[:,n], linewidth = lwidth, color = c1, linestyle = l1)
        dx = 1.8*max(data_315[:,n])
        axes = plt.gca()
        axes.set_ylim([-1.2*max(data_315[:,n]),(len(station_135)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
    else: 
        plt.plot(data_315[:,0], data_315[:,n] + i*dx, linewidth = lwidth, color = c1, linestyle = l1)
    ticks [i] = dx*i

plt.legend(('NE', 'NW', 'SW', 'SE'), prop={'size': 8}, loc='upper left')
axes.set_yticks(ticks)
axes.set_yticklabels(labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
axes.set_ylabel('Distance from sourcs (km)')

# plt.savefig('test.png', dpi = 600)

# %%

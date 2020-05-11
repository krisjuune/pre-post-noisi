# %%
from netCDF4 import Dataset
from obspy.core import Stats, Trace, UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

# %% plot amplitudes - could be useful for thesis! 

# latlon = np.load('/Users/kuangdai/Downloads/sourcegrid.npy')
# plt.figure(dpi=200)
# for i in np.arange(2005):
#     st = 'II.ST' + str(i) + '.RTZ.strain'
#     data = nc.variables[st][:,:]
#     chi = data[:, 0]
#     m = np.max(np.abs(chi))/1e6
#     if (i==22):
#         plt.scatter(latlon[0, i], latlon[1, i], s=5, color=cm.gray(m))
#         print(m)
#         print(chi)
#     else:
#         plt.scatter(latlon[0, i], latlon[1, i], s=5, color=cm.gray(m))
# plt.show()    

# %% get netcdf file
 
 nc = Dataset('/Users/kjoon/Documents/runs/test/output/stations/axisem3d_synthetics.nc')

# %% choose component to plot

 # com='ENZ'
com='RTZ.strain'
data0 = nc.variables['II.ST1.' + com][:,:]
data1 = nc.variables['II.ST2.' + com][:,:]
data2 = nc.variables['II.ST3.' + com][:,:]
data3 = nc.variables['II.ST4.' + com][:,:]

# %% get the potential chi

chi0 = data0[:, 0]
chi1 = data1[:, 0]
chi2 = data2[:, 0]
chi3 = data3[:, 0]

# %% get time stats to trace

time = nc.variables['time_points'][:]
stats = Stats()
stats.starttime = UTCDateTime(time[0])
stats.delta = time[1] - time[0]
stats.npts = len(time)

# %% get data to traces
trace0 = Trace(chi0, header=stats).differentiate().differentiate()/1e+20
trace1 = Trace(chi1, header=stats).differentiate().differentiate()/1e+20
trace2 = Trace(chi2, header=stats).differentiate().differentiate()/1e+20
trace3 = Trace(chi3, header=stats).differentiate().differentiate()/1e+20

# %% 

### plot as sublots, this does mot read the time data... ###
fig, axs = plt.subplots(4, 1)

axs[0].plot(trace0)
axs[1].plot(trace1)
axs[2].plot(trace2)
axs[3].plot(trace3)

plt.tight_layout()

### plot one at a time as an obspy trace ###
# trace0.plot()
# trace1.plot()
# trace2.plot()
# trace3.plot()

# %% get presure
stress = trace.differentiate().differentiate()
stress.plot()
stress.filter('bandpass', freqmin=1/100, freqmax=1/20)
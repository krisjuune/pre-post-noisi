# %% 
from netCDF4 import Dataset
from obspy.core import Stats, Trace, UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

# %% plot sourcegrid and OBS looked at

lonlat = np.load('noisi_stuff/sourcegrid.npy')


# %% plot amplitudes

# idea for later is to have an OBS somewhere in SE and compare with and without bathymetry amplitude snapshots

lonlat = np.load('noisi_stuff/sourcegrid.npy')
run = '362_148'
nc = Dataset('../../../Desktop/flat_Greens/' + run + '/output/stations/axisem3d_synthetics.nc')
h5 = '../noisi/axisem/greens' + 'axisem/greens/II.pres.RTZ.362_148.h5'


plt.figure(dpi=600)
for i in range(len(lonlat)):
    st = 'II.ST' + str(i) + '.RTZ.strain'
    data = nc.variables[st][:,:]
    chi = data[:, 0]
    m = np.max(np.abs(chi))/1e8
    if (i==22):
        plt.scatter(lonlat[0, i], lonlat[1, i], color=cm.gray(m))
        print(m)
        print(chi)
    else:
        plt.scatter(lonlat[0, i], lonlat[1, i], color=cm.gray(m))
plt.show()  

# %%

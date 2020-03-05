from obspy import *
from obspy.clients.fdsn import Client
import numpy as np
import matplotlib.pylab as plt
plt.style.use('ggplot')
from pathlib import Path

# TODO problem with filtering, it is effectively doing nothing atmm 
# perhaps something is wrong with plotting and I am instead plotting 
# the original thingy

# plotting the unfiltered trace works though, so turning it into a
# a trace is fine

# TODO eliminate the grid from the plots completely or add it back 
# onto subplots 1 and 2

# %% Data retrieval 
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

path1 = 'benchmark/processing/raw_data/curvature/geographic/stations'
station = 'ST20' # furthest station in SE
data1 = station_data(path1, station)
data2 = data1.copy()
data3 = data1.copy()

path1_bathy = 'benchmark/processing/raw_data/curvature/geographic_bathymetry/stations' 
station_bathy = 'ST20' # furthest station in SE
data1_bathy = station_data(path1_bathy, station_bathy)
data2_bathy = data1_bathy.copy()
data3_bathy = data1_bathy.copy()

# %% Filtering 
# obtain traces
tr1 = Trace(data1[:,n])
tr1_bathy = Trace(data1_bathy[:,n])
tr2 = Trace(data1[:,n])
tr2_bathy = Trace(data2_bathy[:,n])
tr3 = Trace(data1[:,n])
tr3_bathy = Trace(data3_bathy[:,n])

# set filter parameters
f10 = 0.105                                 # freq for 10s filter
f20 = 0.055                                # freq for 20s filter
corners = 10                              # order of filter
npts = tr2.stats.npts                     # number of samples in the trace
dt = tr2.stats.delta                      # sample interval
fNy = 1. / (2. * dt)                      # Nyquist frequency
time = np.arange(0, npts) * dt            # time axis for plotting
freq = np.linspace(0, fNy, npts // 2 + 1) # frequency axis for plotting

# set desired compenent, 1 radial, 2 transverse, 3 vertical
n = 2

# TODO cool things happen with transverse, bathymetry really has an amplifying effect
# choose vertical for 'regular' plots

# filter woopwoop
tr2.filter('lowpass', freq = f10, corners = corners)
tr2_bathy.filter('lowpass', freq = f10, corners = corners)
tr3.filter('lowpass', freq = f20, corners = corners)
tr3_bathy.filter('lowpass', freq = f20, corners = corners)

# amplitude spectra, BUG - all the amplitude is at freq 0 atm
sp2 = np.fft.rfft(tr2.data)
sp2_bathy = np.fft.rfft(tr2_bathy.data)
sp3 = np.fft.rfft(tr3.data)
sp3_bathy = np.fft.rfft(tr3_bathy.data)

# filter functions
LP10 = 1 / ( 1 + (freq / f10) ** (2 * corners))
LP20 = 1 / ( 1 + (freq / f20) ** (2 * corners))

# %% Plot seismograms
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

plt.figure(1)
# set desired range of data (time plotted)
m = np.arange(0,4000)
lwidth = '0.75'

fig = plt.subplot(311)
plt.plot(data1[:,0], tr1.data, color = 'orange', \
    linewidth = lwidth)
plt.plot(data1_bathy[:,0], tr1_bathy.data, color = 'darkblue', \
    linewidth = lwidth, alpha = 0.70)
plt.legend(('Flat', 'Bathymetry'), \
    prop={'size': 6}, loc='upper right')
plt.title('Unfiltered', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)

fig1 = plt.subplot(312)
plt.plot(data2[:,0], tr2.data, color = 'orange', \
    linewidth = lwidth)
plt.plot(data2_bathy[:,0], tr2_bathy.data, color = 'darkblue', \
    linewidth = lwidth, alpha = 0.70)
plt.legend(('Flat', 'Bathymetry'), \
    prop={'size': 6}, loc='upper right')
plt.title('Filtered at 10s', fontsize = 10)
fig1.axes.get_xaxis().set_visible(False)

plt.subplot(313)
plt.plot(data3[:,0], tr3.data, color = 'orange', \
    linewidth = lwidth)
plt.plot(data3_bathy[:,0], tr3_bathy.data, color = 'darkblue', \
    linewidth = lwidth, alpha = 0.70)
plt.legend(('Flat', 'Bathymetry'), \
    prop={'size': 6}, loc='upper right')
plt.title('Filtered at 20s', fontsize = 10)

plt.subplots_adjust(hspace=0.4)
plt.show()
# plt.savefig('freq_bathy_effect.png', dpi = 600)

# %% Plot amplitude spectra
fx2 = 0.15

plt.plot(freq, LP10, '#ff7f04', linewidth=1.5)
plt.plot(freq, LP20, '#1f77b4', linewidth=1.5)
plt.xlim(0,fx2)
plt.ylim(-0.1,1.1)
plt.title('Filter function', fontsize = 10)
plt.legend(('10s', '20s'), \
    prop={'size': 10}, loc='upper right')
plt.ylabel('amplitude (%)')
plt.xlabel('frequency (Hz)')

plt.savefig('filter_functions.png', dpi = 600)

# Plot saved under thesis plots with corners = 10


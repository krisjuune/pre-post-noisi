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

path1 = '../../runs/benchmark/test_curvature/' + 'Cartesian_3layers' + '/output/stations/'
path2 = '../../runs/benchmark/test_curvature/' + 'Spherical_3layers' + '/output/stations/'
path3 = '../../runs/benchmark/test_curvature/' + 'Geographic_3layers' + '/output/stations/'

station_mid = 'ST4'
data1_mid, time1 = synthetics_data(path1, station_mid)
data2_mid, time2 = synthetics_data(path2, station_mid)
data3_mid, time3 = synthetics_data(path3, station_mid)

station_far = 'ST2'
data1_far, time1 = synthetics_data(path1, station_far)
data2_far, time2 = synthetics_data(path2, station_far)
data3_far, time3 = synthetics_data(path3, station_far)

station_src = 'ST1'
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
plt.plot(time1[:-2570], data1_src[:-2570,n]*(-1)/1e13, color = 'orange', linestyle = '--', linewidth = lwidth)
plt.plot(time2[:-2570], data2_src[:-2570,n]/1e13, color = 'orange', linewidth = lwidth)
plt.plot(time3[:-2570], data3_src[:-2570,n]/1e13, color = '#440154FF', linewidth = lwidth, alpha = 0.95)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('Below source', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)

fig1 = plt.subplot(312)
plt.plot(time1[:-2570], data1_far[:-2570,n]*(-1)/1e13, color = 'orange', linestyle = '--', linewidth = lwidth)
plt.plot(time2[:-2570], data2_far[:-2570,n]/1e13, color = 'orange', linewidth = lwidth)
plt.plot(time3[:-2570], data3_far[:-2570,n]/1e13, color = '#440154FF', linewidth = lwidth, alpha = 0.95)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('85 km', fontsize = 10)
fig1.axes.get_xaxis().set_visible(False)

plt.subplot(313)
plt.plot(time1[:-2570], data1_mid[:-2570,n]*(-1)/1e13, color = 'orange', linestyle = '--', linewidth = lwidth)
plt.plot(time2[:-2570], data2_mid[:-2570,n]/1e13, color = 'orange', linewidth = lwidth)
plt.plot(time3[:-2570], data3_mid[:-2570,n]/1e13, color = '#440154FF', linewidth = lwidth, alpha = 0.95)
plt.legend(('Cartesian', 'Spherical', 'Geographic'), \
    prop={'size': 6}, loc='upper right')
plt.title('170 km', fontsize = 10)
plt

plt.subplots_adjust(hspace=0.4)
plt.savefig('which_curvature_p1e7.png', dpi = 600)
plt.show()

# %% Plot reg vs oversampled TODO look completely different, should be the same
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.filter import *

path4 = '../../runs/benchmark/test_curvature/' + 'Geographic_2.5mesh' + '/output/stations/'
data4_mid, time4 = synthetics_data(path4, station_mid)
data4_far, time4 = synthetics_data(path4, station_far)
data4_src, time4 = synthetics_data(path4, station_src)

path5 = '../../runs/benchmark/test_curvature/' + 'Geographic_doublethick' + '/output/stations/'
data5_mid, time5 = synthetics_data(path5, station_mid)
data5_far, time5 = synthetics_data(path5, station_far)
data5_src, time5 = synthetics_data(path5, station_src)

plt.figure(1)
n = 2
lwidth = '0.5'

fig = plt.subplot(311)
plt.plot(time5, data5_mid[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(time4, data4_mid[:,n], color = 'orange', linewidth = lwidth)
# plt.plot(time3, data3_mid[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('2.5s mesh', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('170 km', fontsize = 10)
fig.axes.get_xaxis().set_visible(False)

fig1 = plt.subplot(312)
plt.plot(time5, data5_far[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(time4, data4_far[:,n], color = 'orange', linewidth = lwidth)
# plt.plot(time3, data3_far[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('2.5s mesh', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('60 km', fontsize = 10)
fig1.axes.get_xaxis().set_visible(False)

plt.subplot(313)
plt.plot(time5, data5_src[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(time4, data4_src[:,n], color = 'orange', linewidth = lwidth)
# plt.plot(time3, data3_src[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('2.5s mesh', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('below source', fontsize = 10)
plt

plt.subplots_adjust(hspace=0.4)
plt.savefig('oversampled.png', dpi = 600)
plt.show()


# %% plot freq spectra for oversampled

duration = 300
Fs4 = 1./(np.mean(time4[1:] - time4[:-1]))
Fs3 = 1./(np.mean(time3[1:] - time3[:-1]))

print(Fs3, Fs4)

# faxis4 = np.fft.rfftfreq(n=duration, d=1./Fs4)
# faxis3 = np.fft.rfftfreq(n=duration, d=1./Fs3)
freq = np.linspace(0, fNy, len(y_f))

spec4 = np.fft.rfft(data4_src, duration)
spec3 = np.fft.rfft(data3_src, duration)

plt.figure(1)
n = 2
lwidth = '0.5'

plt.subplot(211)
plt.plot(faxis4, data4_src[:,n], color = 'darkblue', linewidth = lwidth)
plt.plot(faxis3, data3_src[:,n], color = 'orange', linewidth = lwidth)
plt.legend(('2.5s mesh', '5s mesh'), prop={'size': 6}, loc='upper right')
plt.title('170 km', fontsize = 10)



# %% plot bathy vs flat at 5s, 10s, 20s, 40s
from obspy.signal.filter import * 

path_flat = ['../../runs/benchmark/test_curvature/' + 'Geographic_3layers' + '/output/stations/',
             '../../runs/benchmark/test_curvature/' + 'Geographic_10s' + '/output/stations/',
             '../../runs/benchmark/test_curvature/' + 'Geographic_20s' + '/output/stations/', 
             '../../runs/benchmark/test_curvature/' + 'Geographic_40s' + '/output/stations/']
path_bath = ['../../runs/benchmark/test_curvature/' + 'Geographic_3layers_bathy' + '/output/stations/', 
             '../../runs/benchmark/test_curvature/' + 'Geographic_3layers_bathy_10s' + '/output/stations/', 
             '../../runs/benchmark/test_curvature/' + 'Geographic_3layers_bathy_20s' + '/output/stations/', 
             '../../runs/benchmark/test_curvature/' + 'Geographic_3layers_bathy_40s' + '/output/stations/']

labels_NE = ['6', '7', '8', '9', '10'] # actually NW
labels_SW = ['16', '17', '18', '19', '20'] # actually SE
freqs = ['5', '10', '20', '40']
freqs_labels = ['5 s', '10 s', '20 s', '40 s']
ticks = np.zeros(len(freqs))
c1 = 'dimgrey'
c2 = '#440154FF'
cflat = 'dimgrey'
ls_flat = 'solid'
ls_bath = 'solid'
lw = 0.75
spacing = 2
comp = 2 #ENZ (or RTZ?)
st = 3 #which station in labels

# data_flat_SW, time_SW_flat = synthetics_data(path_flat, 'ST' + labels_SW[2])
# fig = plt.figure()
# plt.plot(time_SW_flat, data_flat_SW[:,2]/1e13)

fig = plt.figure()
for i in range(len(freqs)): 
    station_SW = 'ST' + labels_SW[st]

    if i == 0: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,comp], color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,comp], color=c2, linewidth=lw, linestyle=ls_bath)
        dx = max(data_flat_SW[:,comp])*spacing
        axes = plt.gca()
        axes.set_ylim([-1.5*max(data_flat_SW[:,comp]),(len(freqs)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
        axes.set_xlim([-10,302])

    elif i == 1: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,2]*5+dx*i, color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,2]*5.5+dx*int(i), color=c2, linewidth=lw, linestyle=ls_bath)

    elif i == 2: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,2]*40+dx*i, color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,2]*60+dx*int(i), color=c2, linewidth=lw, linestyle=ls_bath)

    elif i == 3: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,2]*4000+dx*i, color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,2]*5200+dx*int(i), color=c2, linewidth=lw, linestyle=ls_bath)


    ticks[i] = dx*i

plt.legend(('Flat sea-floor', 'Undulating sea-floor'), prop={'size': 8}, bbox_to_anchor=(1.35, 0.57), frameon=False)
axes.set_yticks(ticks)
axes.set_yticklabels(freqs_labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
plt.title('Station 19, SE quadrant')
plt.savefig('freq_bathy_SE.png', dpi = 600)
plt.show()




fig = plt.figure()
for i in range(len(freqs)): 
    station_SW = 'ST' + labels_NE[st]

    if i == 0: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,comp], color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,comp], color=c2, linewidth=lw, linestyle=ls_bath)
        dx = max(data_flat_SW[:,comp])*spacing
        axes = plt.gca()
        axes.set_ylim([-1.5*max(data_flat_SW[:,comp]),(len(freqs)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
        axes.set_xlim([-10,302])

    elif i == 1: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,2]*5+dx*i, color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,2]*5.5+dx*int(i), color=c2, linewidth=lw, linestyle=ls_bath)

    elif i == 2: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,2]*40+dx*i, color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,2]*60+dx*int(i), color=c2, linewidth=lw, linestyle=ls_bath)

    elif i == 3: 
        data_flat_SW, time_SW_flat = synthetics_data(path_flat[i], station_SW)
        data_bath_SW, time_SW_bath = synthetics_data(path_bath[i], station_SW)
        # data_flat_SW = lowpass(data_flat_SW, 1/int(freqs[i]), 50, corners=9)
        plt.plot(time_SW_flat, data_flat_SW[:,2]*4000+dx*i, color=cflat, linewidth=lw, linestyle=ls_flat)
        plt.plot(time_SW_bath, data_bath_SW[:,2]*5500+dx*int(i), color=c2, linewidth=lw, linestyle=ls_bath)


    ticks[i] = dx*i

plt.legend(('Flat sea-floor', 'Undulating sea-floor'), prop={'size': 8}, bbox_to_anchor=(1.35, 0.57), frameon=False)
axes.set_yticks(ticks)
axes.set_yticklabels(freqs_labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
plt.title('Station 9, NW quadrant')
plt.savefig('freq_bathy_NW.png', dpi = 600)
plt.show()


# %% Love wave test nono - needs noisi of course.... 
# path1 = '../../runs/Love_wave/' + 'flat' + '362_160/output/stations/'
# path2 = '../../runs/Love_wave/' + 'bathy' + '362_162/output/stations/'

# station_mid = 'ST500'
# data1_mid, time1 = synthetics_data(path1, station_mid)
# data2_mid, time2 = synthetics_data(path2, station_mid)
# data3_mid, time3 = synthetics_data(path3, station_mid)


# %%
import numpy as np 
import obspy # allows to filter as lowpass
from pathlib import Path
# from benchmark.functions import station_data
import matplotlib.pyplot as plt 
    

station_NW = ('ST6', 'ST7', 'ST8', 'ST9', 'ST10')
station_SE = ('ST16', 'ST17', 'ST18', 'ST19', 'ST20')
ticks = np.zeros(len(station_NW))
labels = ['40', '85', '130', '170', '210']



path_flat = '../../runs/benchmark/test_curvature/' + 'Geographic_3layers' + '/output/stations/'
path_bath = '../../runs/benchmark/test_curvature/' + 'Geographic_3layers_bathy' + '/output/stations/'
file_flat = 'NW_moveout.png'
file_bath = 'SE_moveout.png'


fig1 = plt.figure()
for i in np.arange(len(ticks)):
    data_NW, time_NW = synthetics_data(path_flat, station_NW [i])
    data_SE, time_SE = synthetics_data(path_bath, station_NW [i])

    # choose component to be plotted
    n = 2
    lwidth = 0.75
    c1 = 'dimgrey'
    c2 = '#440154FF'
    l1 = 'solid'
    l2 = 'solid'

    if i == 0: 
        plt.plot(time_NW, data_NW[:,n]/data_NW[:,n].max(), linewidth = lwidth, color = c1, linestyle = l1, alpha = 0.65)
        plt.plot(time_SE, data_SE[:,n]/data_SE[:,n].max(), linewidth = lwidth, color = c2, linestyle = l2, alpha = 0.85)
        dx = 2.0*max(data_NW[:,n]/data_NW[:,n].max())
        axes = plt.gca()
        axes.set_ylim([-1.2*max(data_NW[:,n]/data_NW[:,n].max()),(len(station_SE)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
    else: 
        plt.plot(time_NW, data_NW[:,n]/data_NW[:,n].max() + i*dx, linewidth = lwidth, color = c1, linestyle = l1, alpha = 0.65)
        plt.plot(time_SE, data_SE[:,n]/data_SE[:,n].max() + i*dx, linewidth = lwidth, color = c2, linestyle = l2, alpha = 0.85)
    ticks [i] = dx*i

plt.legend(('flat', 'undulating'), prop={'size': 8}, loc='upper left')
axes.set_yticks(ticks)
axes.set_yticklabels(labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
axes.set_ylabel('Distance from source (km)')
plt.title('NW quadrant')

plt.savefig(file_flat, dpi = 600)
plt.show


fig2 = plt.figure()
for i in np.arange(len(ticks)):
    data_NW, time_NW = synthetics_data(path_flat, station_SE [i])
    data_SE, time_SE = synthetics_data(path_bath, station_SE [i])

    # choose component to be plotted
    n = 2
    lwidth = 0.75
    c1 = 'dimgrey'
    c2 = '#440154FF'
    l1 = 'solid'
    l2 = 'solid'

    if i == 0: 
        plt.plot(time_NW, data_NW[:,n]/data_NW[:,n].max(), linewidth = lwidth, color = c1, linestyle = l1, alpha = 0.65)
        plt.plot(time_SE, data_SE[:,n]/data_SE[:,n].max(), linewidth = lwidth, color = c2, linestyle = l2, alpha = 0.85)
        dx = 2.0*max(data_NW[:,n]/data_NW[:,n].max())
        axes = plt.gca()
        axes.set_ylim([-1.2*max(data_NW[:,n]/data_NW[:,n].max()),(len(station_SE)-0.3)*dx])
        axes.get_yaxis().set_visible(False)
    else: 
        plt.plot(time_NW, data_NW[:,n]/data_NW[:,n].max() + i*dx, linewidth = lwidth, color = c1, linestyle = l1, alpha = 0.65)
        plt.plot(time_SE, data_SE[:,n]/data_SE[:,n].max() + i*dx, linewidth = lwidth, color = c2, linestyle = l2, alpha = 0.85)
    ticks [i] = dx*i

plt.legend(('flat', 'undulating'), prop={'size': 8}, loc='upper left')
axes.set_yticks(ticks)
axes.set_yticklabels(labels)
axes.get_yaxis().set_visible(True)
axes.set_xlabel('Time (s)')
axes.set_ylabel('Distance from source (km)')
plt.title('SE quadrant')

plt.savefig(file_bath, dpi = 600)
plt.show





# %%

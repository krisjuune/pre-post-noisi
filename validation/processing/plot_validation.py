import numpy as np 
import obspy as obs 
from pathlib import Path
import matplotlib.pyplot as plt

# %% Get the data
def station_data(path, station):
    """
    Function that retrieves the seismic data from station
    'station' given the relative path 'path', both inputs 
    are strings. This works for II type (not IU) stations.
    Function returns the data array. 
    """
    from pathlib import Path
    import numpy as np
    path = Path(path)
    file2open = 'II.' + station + '.RTZ.ascii'
    # file handle
    file2open = path/file2open 
    # Open file and retrieve data
    data = open(file2open, 'r')
    data = data.read()
    data = np.fromstring(data, dtype = float, sep=' ')
    m = len(data)/4
    data = data.reshape(m,4)
    return(data)

path1 = 'validation/processing/raw_data/sim1/'
path2 = 'validation/processing/raw_data/sim2/'
path3 = 'validation/processing/raw_data/sim3/'
station = 'LAT5'
data1 = station_data(path1, station)
data2 = station_data(path2, station)
data3 = station_data(path3, station)

# %% Plot seismograms
cmap = 'viridis'
plt.plot(data1[0:400,0], data1[0:400,4], cmap=cmap)
import numpy as np 
from pathlib import Path
import matplotlib.pyplot as plt
# from benchmark.processing.processing_functions import \
#     station_data

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

stations_lat = ['LAT2', 'LAT4', 'LAT6', 'LAT8', 'LAT10']
stations_lon = ['LON2', 'LON4', 'LON6', 'LON8', 'LON10']
path = 'benchmark/processing/raw_data/sim1'

suva_data = station_data(path, stations_lat [1])

def plot_dist_time(path, stations, save = False, \
    *filename, vertical = True, transverse = False, \
    radial = False, data_points = 500):
    """
    Function to plot seismograms stack atop one another as
    distance vs time, using the function station_data. 
    Inputs are a string array with station names to be 
    plotted, a string with the path to the station files, 
    and a boolean statement whether the plot is to be saved. 
    If set save to True, add the desired filename as a str. 
    Choose the component to be plotted, either vertical 
    or transverse. Vertical chosen by default. 
    """
    fig = plt.figure(figsize = (7,7), dpi = 1200)

    for i in np.arange(len(stations)):
        # Get the data from all the stations
        data = station_data(path, stations[i])
        time = data[0:data_points,0]
        # Retrieve the desired comoponent data
        if vertical == True:
            displacement = data[0:data_points,3]
        elif transverse == True: 
            displacement = data[0:data_points,2]
        elif radial == True:
            displacement = data[0:data_points,1]
        else: 
            print('Error, choose displacement component')
        
        # Begin plotting
        if i == 0:
            plt.plot(time, displacement, linewidth = 1, \
                linestyle = '-', color = 'orange')
            dx = 0.8*max(displacement)
    return fig

plot_dist_time(path, stations_lat)

    


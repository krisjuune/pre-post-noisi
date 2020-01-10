import numpy as np 
from pathlib import Path
import matplotlib.pyplot as plt
from benchmark.processing.plot_validation import station_data

stations_lat = ['LAT2', 'LAT4', 'LAT6', 'LAT8', 'LAT10']
stations_lon = ['LON2', 'LON4', 'LON6', 'LON8', 'LON10']
path = benchmark/processing/raw_data/sim1

def plot_dist_time(path, stations, save = False, \
    *filename, vertical = True, transverse = False, \
    radial = False):
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

    for i in np.arange(len(stations)):
        
        # Get the data from all the stations
        data = station_data(path, stations[i])
        time = data[:,0]
        # Retrieve the desired comoponent data
        if vertical == True:
            displacement = data[:,3]
        elif transverse == True: 
            displacement = data[:,2]
        elif radial == True:
            displacement = data[:,1]
        else: 
            print('Error, choose displacement component')
        
        # Begin plotting
        if i == 0:
            plt(time, displacement, linewidth = 1, linestyle = '-', \
                color = 'orange')
            dx = 0.8*max(displacement)
        plt.show()

    


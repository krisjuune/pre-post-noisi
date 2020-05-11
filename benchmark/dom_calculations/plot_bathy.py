# need matplotlib, basemap, anything else?

# %% Truncate bathymetry data
import pickle 
import numpy as np
import numpy.ma as ma
from pathlib import Path
from coordinate_transformation.functions.get_domain import \
    find_nearest, truncate_domain

# Load the dumped variables
path = 'coordinate_transformation/variables/'
path = Path(path)

with open(path / 'lat_Prt', 'rb') as f:
    lat_Prt = pickle.load(f)
lat_Prt = np.ma.getdata(lat_Prt)
with open(path / 'lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)
with open(path / 'bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)

# Define domain (within 35.5...39.5N, -14...-19W)
lat_max = 39.5
lat_min = 35.5
lon_max = -14
lon_min = -19
bounds = [lat_max, lat_min, lon_max, lon_min]

(lat_dom, lon_dom, elevation_dom) = truncate_domain(lat_Prt, \
    lon_Prt, bathy_Prt, bounds)

elevation_dom = elevation_dom/1000 #plot elevation in km

# %% Plot bathymetry using Basemap, info from ObsPy basemap plots
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.basemap import Basemap

def plot_geographic(lat, lon, data, filename, \
    lat_max = 39.5, lat_min = 35.5, lon_max = -14, \
    lon_min = -19, cbar_label = 'Bathymetry (km)'):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig = Basemap(projection = 'mill', llcrnrlat = lat_min, \
        urcrnrlat = lat_max, llcrnrlon = lon_min, \
        urcrnrlon = lon_max, resolution = 'c')
    fig.drawmapboundary()
    # Draw a lon/lat grid (20 lines for an interval of one degree)
    fig.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
        labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 2])
    fig.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
        labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 2])

    # Add elevation data to map
    cmap = 'viridis'
    Lon, Lat = np.meshgrid(lon, lat)
    fig.pcolormesh(Lon, Lat, data, latlon = True, \
        cmap = cmap)
    # Colorbar construction
    i = ax.imshow(data, interpolation='nearest')
    cbar = fig.colorbar(i, shrink = 0.5, aspect = 5)
    cbar.set_label(cbar_label, rotation = 270, labelpad=15, y=0.45)

    plt.savefig(filename, dpi = 600)
    plt.show()

plot_geographic(lat_dom, lon_dom, elevation_dom, 'test.png')

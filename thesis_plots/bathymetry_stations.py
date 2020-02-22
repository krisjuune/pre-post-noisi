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

# %% Plot bathymetry using Basemap, Miller projection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.basemap import Basemap
from coordinate_transformation.functions.transform import \
    cartesian_to_geographic

filename = 'stations.png'
lat_max = 39.5 
lat_min = 35.5 
lon_max = -14
lon_min = -19 
cbar_label = 'Bathymetry (km)'

##### PLOT BATHYMETRY #####
map = plt.figure()
ax = map.add_subplot(1, 1, 1)
map = Basemap(projection = 'mill', llcrnrlat = lat_min, \
    urcrnrlat = lat_max, llcrnrlon = lon_min, \
    urcrnrlon = lon_max, resolution = 'c')
map.drawmapboundary()
# Draw a lon/lat grid (20 lines for an interval of one degree)
map.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
    labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 2], fontsize = 8)
map.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
    labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 2], fontsize = 8)

# Add elevation data to map
cmap = 'viridis'
Lon, Lat = np.meshgrid(lon_dom, lat_dom)
map.pcolormesh(Lon, Lat, elevation_dom, latlon = True, \
    cmap = cmap)
# Colorbar construction, TODO set cbar fontsize to 8 
i = ax.imshow(elevation_dom, interpolation='nearest')
cbar = map.colorbar(i, shrink = 0.5, aspect = 5)
cbar.set_label(cbar_label, rotation = 270, labelpad=15, y=0.45, \
    fontsize = 8)

##### ADD STATIONS TO PLOT #####
x_stations = np.array([0, \
    28.3, 60.1, 91.9, 120.2, 149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    28.3, 60.1, 91.9, 120.2, 149.5]) * 1000
y_stations = np.array([0, \
    28.3, 60.1, 91.9, 120.2, 149.5, \
    28.3, 60.1, 91.9, 120.2, 149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5]) * 1000

# get latitude and longitude of stations for plotting
lat_stations, lon_stations = cartesian_to_geographic(x_stations, \
    y_stations)
name = np.zeros(len(x_stations), int)
xpt = np.zeros(len(x_stations))
ypt = np.zeros(len(x_stations))
x_dist = 8000
y_dist = -6000

# add stations to plot, looping over each station
for i in np.arange(len(x_stations)): 
    map.plot(lon_stations[i], lat_stations[i], 'bv', \
        latlon = True, markersize = 4)
    name[i] = str(i) # station number
    xpt[i], ypt[i] = map(lon_stations[i], lat_stations[i])
    plt.text(xpt[i]+x_dist, ypt[i]+y_dist, name[i], fontsize = 8, \
        color = 'darkblue')

plt.savefig(filename, dpi = 600)
plt.show()

# %%

# %%

import pickle 
import numpy as np
import numpy.ma as ma
from pathlib import Path
from coordinate_transformation.functions.get_domain import \
    find_nearest, truncate_domain
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap
from coordinate_transformation.functions.transform import \
    cartesian_to_geographic

# %% Truncate bathymetry data

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
map.drawmapboundary(linewidth = 0.8)
# Draw a lon/lat grid (20 lines for an interval of one degree)
map.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
    labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
map.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
    labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')

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
    map.plot(lon_stations[i], lat_stations[i], color = 'orange', marker = 'v', \
        latlon = True, markersize = 4)
    name[i] = str(i) # station number
    xpt[i], ypt[i] = map(lon_stations[i], lat_stations[i])
    plt.text(xpt[i]+x_dist, ypt[i]+y_dist, name[i], fontsize = 8, \
        color = '#ff7f04')

plt.savefig(filename, dpi = 600)
plt.show()

# %% plot noisi NS test

color_bg = [62/255, 74/255, 137/255]

filename = 'noisi_NS_test_setup.png'
lat_max = 39.5 
lat_min = 35.5 
lon_max = -14
lon_min = -19 
cbar_label = 'Bathymetry (km)'
cmap = 'viridis'

##### PLOT BACKGROUND #####
map = plt.figure()
ax = map.add_subplot(1, 1, 1)
map = Basemap(projection = 'mill', llcrnrlat = lat_min, \
    urcrnrlat = lat_max, llcrnrlon = lon_min, \
    urcrnrlon = lon_max, resolution = 'c')
map.drawmapboundary(fill_color=color_bg, linewidth = 0.8)
# Draw a lon/lat grid (20 lines for an interval of one degree)
map.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
    labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
map.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
    labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
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
lat_stations, lon_stations = cartesian_to_geographic(x_stations, y_stations)
name = np.zeros(len(x_stations), int)
xpt = np.zeros(len(x_stations))
ypt = np.zeros(len(x_stations))
x_dist = 8000
y_dist = -6000

# add stations to plot, looping over each station
for i in np.arange(len(x_stations)): 
    map.plot(lon_stations[i], lat_stations[i], color = 'k', marker = 'v', \
        latlon = True, markersize = 4, alpha = 0.8)
    name[i] = str(i) # station number
    xpt[i], ypt[i] = map(lon_stations[i], lat_stations[i])
    plt.text(xpt[i]+x_dist, ypt[i]+y_dist, name[i], fontsize = 8, \
        color = 'k')


##### ADD CIRCLES TO PLOT #####
lon_circlesN = np.array([lon_stations[16], lon_stations[13]])
lat_circlesN = np.array([lat_stations[16], lat_stations[13]])
lon_circlesS = np.array([lon_stations[6], lon_stations[3]])
lat_circlesS = np.array([lat_stations[6], lat_stations[3]])
nameN = ['16', '13']
nameS = ['6', '3']
xptN = np.zeros(len(lon_circlesN))
yptN = np.zeros(len(lon_circlesN))
xptS = np.zeros(len(lon_circlesS))
yptS = np.zeros(len(lon_circlesS))
x_dist = 8000
y_dist = -6000
for i in np.arange(len(lat_circlesN)):
    map.plot(lon_circlesN[i], lat_circlesN[i], color = 'yellowgreen', marker = 'v', \
        latlon = True, markersize = 4, alpha = 1)
    xptN[i], yptN[i] = map(lon_circlesN[i], lat_circlesN[i])
    plt.text(xptN[i]+x_dist, yptN[i]+y_dist, nameN[i], fontsize = 8, \
        color = 'limegreen')
for i in np.arange(len(lat_circlesS)):
    map.plot(lon_circlesS[i], lat_circlesS[i], color = 'orange', marker = 'v', \
        latlon = True, markersize = 4, alpha = 1)
    xptS[i], yptS[i] = map(lon_circlesS[i], lat_circlesS[i])
    plt.text(xptS[i]+x_dist, yptS[i]+y_dist, nameS[i], fontsize = 8, \
        color = '#ff7f04')

##### ADD NOISE SOURCES TO PLOT #####
lat_sources = np.array([38.7, 36.29967440692094])
lon_sources = np.array([-15.5, -17.5])
name_sources = ['NE', 'SW']
xpt = np.zeros(len(lat_sources))
ypt = np.zeros(len(lat_sources))
map.plot(lon_sources[0], lat_sources[0], color = 'yellowgreen', marker = 'o', \
    latlon = True, markersize = 17, alpha=0.1)
map.plot(lon_sources[0], lat_sources[0], color = 'yellowgreen', marker = 'o', \
    latlon = True, markersize = 10, alpha=0.5)
map.plot(lon_sources[0], lat_sources[0], color = 'yellowgreen', marker = 'o', \
    latlon = True, markersize = 3, alpha=1)
xpt[0], ypt[0] = map(lon_sources[0], lat_sources[0])
plt.text(xpt[0]+14000, ypt[0]+7000, name_sources[0], fontsize = 12, \
        color = 'limegreen')
map.plot(lon_sources[1], lat_sources[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 17, alpha=0.1)
map.plot(lon_sources[1], lat_sources[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 10, alpha=0.5)
map.plot(lon_sources[1], lat_sources[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 3, alpha=1)
xpt[1], ypt[1] = map(lon_sources[1], lat_sources[1])
plt.text(xpt[1]+14000, ypt[1]-30000, name_sources[1], fontsize = 12, \
        color = '#ff7f04')

plt.savefig(filename, dpi = 600)
plt.show()


# %% plot noisi cnt test setup

color_bg = [62/255, 74/255, 137/255]

filename = 'noisi_cnt_test_setup.png'
lat_max = 39.5 
lat_min = 35.5 
lon_max = -14
lon_min = -19 
cbar_label = 'Bathymetry (km)'
cmap = 'viridis'

##### PLOT BACKGROUND #####
map = plt.figure()
ax = map.add_subplot(1, 1, 1)
map = Basemap(projection = 'mill', llcrnrlat = lat_min, \
    urcrnrlat = lat_max, llcrnrlon = lon_min, \
    urcrnrlon = lon_max, resolution = 'c')
map.drawmapboundary(fill_color=color_bg, linewidth = 0.8)
# Draw a lon/lat grid (20 lines for an interval of one degree)
map.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
    labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
map.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
    labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
# Colorbar construction, TODO set cbar fontsize to 8 
i = ax.imshow(elevation_dom, interpolation='nearest')
cbar = map.colorbar(i, shrink = 0.5, aspect = 5)
cbar.set_label(cbar_label, rotation = 270, labelpad=15, y=0.45, \
    fontsize = 8)





##### ADD STATIONS TO PLOT #####
x_stations = np.array([28.3, 60.1, 91.9, 120.2, 149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    28.3, 60.1, 91.9, 120.2, 149.5]) * 1000
y_stations = np.array([28.3, 60.1, 91.9, 120.2, 149.5, \
    28.3, 60.1, 91.9, 120.2, 149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5]) * 1000

# get latitude and longitude of stations for plotting
lat_stations, lon_stations = cartesian_to_geographic(x_stations, y_stations)
name = np.zeros(len(x_stations), int)
xpt = np.zeros(len(x_stations))
ypt = np.zeros(len(x_stations))
x_dist = 8000
y_dist = -6000

# add stations to plot, looping over each station
for i in np.arange(len(x_stations)): 
    map.plot(lon_stations[i], lat_stations[i], color = 'k', marker = 'v', \
        latlon = True, markersize = 4, alpha = 0.8)
    name[i] = str(i+1) # station number
    xpt[i], ypt[i] = map(lon_stations[i], lat_stations[i])
    plt.text(xpt[i]+x_dist, ypt[i]+y_dist, name[i], fontsize = 8, \
        color = 'k')



##### ADD CIRCLES TO PLOT #####
lon_circlesS = np.array([lon_stations[4], lon_stations[9], lon_stations[14], lon_stations[19]])
lat_circlesS = np.array([lat_stations[4], lat_stations[9], lat_stations[14], lat_stations[19]])
nameS = ['5', '10', '15', '20']
xptS = np.zeros(len(lon_circlesS))
yptS = np.zeros(len(lon_circlesS))
x_dist = 8000
y_dist = -6000
for i in np.arange(len(lat_circlesS)):
    map.plot(lon_circlesS[i], lat_circlesS[i], color = 'orange', marker = 'v', \
        latlon = True, markersize = 4, alpha = 1)
    xptS[i], yptS[i] = map(lon_circlesS[i], lat_circlesS[i])
    plt.text(xptS[i]+x_dist, yptS[i]+y_dist, nameS[i], fontsize = 8, \
        color = '#ff7f04')



##### ADD NOISE SOURCE TO PLOT #####
lat_sources = 37.5
lon_sources = -16.5
map.plot(lon_sources, lat_sources, color = 'orange', marker = 'o', \
    latlon = True, markersize = 17, alpha=0.1)
map.plot(lon_sources, lat_sources, color = 'orange', marker = 'o', \
    latlon = True, markersize = 10, alpha=0.5)
map.plot(lon_sources, lat_sources, color = 'orange', marker = 'o', \
    latlon = True, markersize = 3, alpha=1)

plt.savefig(filename, dpi = 600)
plt.show()



# %% sideways test flat and bathy plots
filename_flat = 'flat_setup.png'
filename_bathy = 'bathy_setup.png'

lat_max = 39.5 
lat_min = 35.5 
lon_max = -14
lon_min = -19 
cbar_label = 'Bathymetry (km)'
source=[-14.80981, 38.84655]
color_bg = [62/255, 74/255, 137/255]

x_stations = np.array([0, \
    -28.3, -60.1, -91.9, -120.2, -149.5, \
    28.3, 60.1, 91.9, 120.2, 149.5]) * 1000
y_stations = np.array([0, \
    28.3, 60.1, 91.9, 120.2, 149.5, \
    -28.3, -60.1, -91.9, -120.2, -149.5]) * 1000

st_name = ['0', '6', '7', '8', '9', '10', '16', '17', '18', '19', '20']

# get latitude and longitude of stations for plotting
lat_stations, lon_stations = cartesian_to_geographic(x_stations, y_stations)
xpt = np.zeros(len(x_stations))
ypt = np.zeros(len(x_stations))
x_dist = 8000
y_dist = -6000


# %%
##### PLOT WITH BATHYMETRY #####
map = plt.figure()
ax = map.add_subplot(1, 1, 1)
map = Basemap(projection = 'mill', llcrnrlat = lat_min, \
    urcrnrlat = lat_max, llcrnrlon = lon_min, \
    urcrnrlon = lon_max, resolution = 'c')
map.drawmapboundary(linewidth = 0.8)
# Draw a lon/lat grid (20 lines for an interval of one degree)
map.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
    labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
map.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
    labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
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
# Add stations
for i in np.arange(len(lon_stations)): 
    map.plot(lon_stations[i], lat_stations[i], color = 'orange', marker = 'v', \
        latlon = True, markersize = 4)
    xpt[i], ypt[i] = map(lon_stations[i], lat_stations[i])
    plt.text(xpt[i]+x_dist, ypt[i]+y_dist, st_name[i], fontsize = 8, \
        color = '#ff7f04')
# Add source
map.plot(source[0], source[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 63, alpha=0.1)
map.plot(source[0], source[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 37, alpha=0.2)
map.plot(source[0], source[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 12, alpha=0.5)
# Save and show
plt.savefig(filename_bathy, dpi = 600)
plt.show()

# %%
##### PLOT WITHOUT BATHYMETRY #####
map = plt.figure()
ax = map.add_subplot(1, 1, 1)
map = Basemap(projection = 'mill', llcrnrlat = lat_min, \
    urcrnrlat = lat_max, llcrnrlon = lon_min, \
    urcrnrlon = lon_max, resolution = 'c')
map.drawmapboundary(fill_color=color_bg, linewidth = 0.8)
# Draw a lon/lat grid (20 lines for an interval of one degree)
map.drawparallels(np.linspace(lat_min, lat_max, num = 5), \
    labels=[1, 0, 0, 0], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
map.drawmeridians(np.arange(round(lon_min), round(lon_max), 1), \
    labels=[0, 0, 0, 1], fmt="%.2f", dashes=[2, 0], fontsize = 8, linewidth = 0.5, color = 'dimgrey')
# Add elevation data to map
cmap = 'viridis'
i = ax.imshow(elevation_dom, interpolation='nearest')
cbar = map.colorbar(i, shrink = 0.5, aspect = 5)
cbar.set_label(cbar_label, rotation = 270, labelpad=15, y=0.45, \
    fontsize = 8)
# Add stations
for i in np.arange(len(lon_stations)): 
    map.plot(lon_stations[i], lat_stations[i], color = 'orange', marker = 'v', \
        latlon = True, markersize = 4)
    xpt[i], ypt[i] = map(lon_stations[i], lat_stations[i])
    plt.text(xpt[i]+x_dist, ypt[i]+y_dist, st_name[i], fontsize = 8, \
        color = '#ff7f04')
# Add source
map.plot(source[0], source[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 63, alpha=0.1)
map.plot(source[0], source[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 37, alpha=0.2)
map.plot(source[0], source[1], color = 'orange', marker = 'o', \
    latlon = True, markersize = 12, alpha=0.5)
# Save and show
plt.savefig(filename_flat, dpi = 600)
plt.show()


# %%

# need matplotlib, basemap, anything else?

#%% Truncate bathymetry data
import pickle 
import numpy as np
import numpy.ma as ma
from transformation_functions.get_domain \
    import find_nearest, truncate_domain

# Load the dumped variables
with open('coordinate_transformation/lat_Prt', 'rb') as f:
    lat_Prt = pickle.load(f)
lat_Prt = np.ma.getdata(lat_Prt)
with open('coordinate_transformation/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)
with open('coordinate_transformation/bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)

# Define domain
lat_max = 39
lat_min = 37
lon_max = -16.5
lon_min = -19.5
bounds = [lat_max, lat_min, lon_max, lon_min]

(lat_dom, lon_dom, elevation_dom) = truncate_domain(lat_Prt, \
    lon_Prt, bathy_Prt, bounds)

elevation_dom = elevation_dom/1000 #plot elevation in km

#%% Plot bathymetry using Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.basemap import Basemap

map = Basemap(projection = 'mill', llcrnrlat = 37, \
    urcrnrlat = 39, llcrnrlon = -16.5, urcrnrlon = -19, 
    resolution = 'c')
map.drawmapboundary()
ax = plt.gca() 
fig = plt.gcf()

# Add elevation data to map
cmap = 'viridis'
Lon_dom, Lat_dom = np.meshgrid(lon_dom, lat_dom)
map.pcolormesh(Lon_dom, Lat_dom, elevation_dom, latlon = True, \
    cmap = cmap)
# Colorbar construction
cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
cb = ColorbarBase(cax,cmap=cmap, orientation='horizontal')
cb.ax.set_xlabel('Bathymetry in chosen domain (km)')

plt.axes(ax)
plt.title('wad is dis')
plt.show()

# %%

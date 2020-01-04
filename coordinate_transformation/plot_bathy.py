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
    lon_Prt, raw_elevation, bounds)

#%% Plot bathymetry using Basemap
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

m = Basemap(projection = 'mill', llcrnrlat = 37, \
    urcrnrlat = 39, llcrnrlon = -16.5, urcrnrlon = -19, 
    resolution = 'c')

m.drawmapboundary()

# Add elevation data to map
Lat_dom, Lon_dom = np.meshgrid(lat_dom, lon_dom)
map.pcolormesh(Lat_dom, Lon_dom, elevation_dom)

plt.title('Chosen domain')
plt.show()

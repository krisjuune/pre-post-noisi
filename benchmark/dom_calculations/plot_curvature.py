# %% Plot curvature using Basemap plots
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.basemap import Basemap

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
cmap = 'inferno'
Lon_dom, Lat_dom = np.meshgrid(lon_dom, lat_dom)
surface_sphere = np.transpose(surface_sphere)
fig.pcolormesh(Lon_dom, Lat_dom, surface_sphere, latlon = True, \
    cmap = cmap)
# Colorbar construction
i = ax.imshow(elevation_dom, interpolation='nearest')
cbar = fig.colorbar(i, shrink = 0.5, aspect = 5)
cbar.set_label('Curvature (km)', rotation = 270, labelpad=15, y=0.45)

plt.show()


# %%

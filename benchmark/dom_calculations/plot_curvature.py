# %% Plot curvature using Basemap plots
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from math import pi
from coordinate_transformation.functions.transform \
    import geographic_to_geocentric, get_cartesian_distance, \
        radius_cnt, wgs84
from coordinate_transformation.functions.domain \
    import get_variable
from benchmark.functions \
    import plot_curvature

surface_sphere = get_variable('surface_sphere', \
    'coordinate_transformation/variables/')
bottom_sphere = get_variable('bottom_sphere', \
    'coordinate_transformation/variables/')
surface_ellipsoid = get_variable('surface_ellipsoid', \
    'coordinate_transformation/variables/')
bottom_ellipsoid = get_variable('bottom_ellipsoid', \
    'coordinate_transformation/variables/')
lon_dom = get_variable('lon_dom', \
    'coordinate_transformation/variables/')
lat_dom = get_variable('lat_dom', \
    'coordinate_transformation/variables/')

# plot_curvature(lat_dom, lon_dom, surface_sphere, \
#     filename = 'surface_sphere_curvature')

plot_curvature(lat_dom, lon_dom, surface_ellipsoid)

plot_curvature(lat_dom, lon_dom, bottom_ellipsoid)

# %% THIS WORKS 
# test plotting outside function
# Transform lat, lon to be centered around the N Pole
src_lat = 37.5
src_lon = -16.5
cbar_label = 'Curvature (km)'
(x, y) = get_cartesian_distance(lon_dom, lat_dom, \
    src_lon, src_lat)
x, y = np.meshgrid(x, y)
x = x.transpose()
y = y.transpose()

# Create figure handle
fig = plt.figure()
ax = plt.gca(projection = '3d')
# TODO how to scale the axes, so z not so exaggerated
# Plot 
surf = ax.plot_surface(x, y, surface_sphere*(-1), \
    cmap = 'viridis')
# Add colorbar
cbar = fig.colorbar(surf, shrink = 0.5, aspect = 5)
cbar.set_label(cbar_label, rotation = 270, labelpad=15, \
    y=0.45)


# %%

# %% Plot curvature using Basemap plots
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from coordinate_transformation.functions.get_rotation \
    import get_cartesian_distance
from coordinate_transformation.functions.get_spherical \
    import geographic_to_geocentric

def plot_curvature(lat, lon, curvature, src_lat = 37.5, \
    src_lon = -16.5, cbar_label = 'Curvature (km)', \
    filename = 'noname'):
    """
    Function to plot a 3d surface once transformed 
    into Cartesian distances. BUG SOMEWHERE, FIX
    Plots fine with transpose outside the function. 
    """
    
    # Transform lat, lon to be centered around the N Pole
    (x, y) = get_cartesian_distance(lon, lat, \
        src_lat, src_lon)
    x, y = np.meshgrid(x, y)
    # x = x.transpose()
    # y = y.transpose()

    # Create figure handle
    fig = plt.figure()
    ax = plt.gca(projection = '3d')
    # TODO how to scale the axes, so z not so exaggerated
    # Plot 
    surf = ax.plot_surface(x, y, curvature*(-1), \
        cmap = 'viridis')
    # Add colorbar
    cbar = fig.colorbar(surf, shrink = 0.5, aspect = 5)
    cbar.set_label(cbar_label, rotation = 270, labelpad=15, \
        y=0.45)

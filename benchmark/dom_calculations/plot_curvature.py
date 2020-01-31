# %% Plot curvature using Basemap plots
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.basemap import Basemap
from benchmark.dom_calculations.functions import \
    plot_geographic

########## does not work ############
############# bollocks ##############

plot_geographic(lat_dom, lon_dom, surface_sphere, \
    'curv_test.png', cbar_label='Curvature (km)')

# %%

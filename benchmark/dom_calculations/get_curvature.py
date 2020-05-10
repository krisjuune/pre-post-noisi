# %% Load variables
import numpy as np 
import numpy.ma as ma
import pickle 
from pathlib import Path
from math import pi
from coordinate_transformation.functions.get_domain \
    import find_nearest, truncate_domain

with open('coordinate_transformation/variables/lon_Prt', 'rb') as f:
    lon_Prt = pickle.load(f)
lon_Prt = np.ma.getdata(lon_Prt)

with open('coordinate_transformation/variables/lat_Prt', 'rb') as f:
    lat_Prt = pickle.load(f)
lat_Prt = np.ma.getdata(lat_Prt)

with open('coordinate_transformation/variables/bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)


# Make the domain below for calculating curvature 
# slightly bigger to have enough data (along lon)
# That way there is enough data here to cover all
lat_max = 40.0
lat_min = 35.0
lon_max = -13.5
lon_min = -19.5
bounds = [lat_max, lat_min, lon_max, lon_min]

(lat_dom, lon_dom, bathy_dom) = truncate_domain(lat_Prt, \
    lon_Prt, bathy_Prt, bounds)

# %% Get arrays with distances to the curve for sphere
import numpy as np 
from math import pi
from benchmark.functions import get_curvature

surface_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6370287.273)
ocean_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6365567.273)
Moho_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6358117.273)
bottom_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6270287.273)

# %% Save curvature values for the spherical case
from netCDF4 import Dataset
import numpy as np 
import datetime as dt
from coordinate_transformation.functions.transform \
    import geographic_to_geocentric, wgs84, radius_cnt
from coordinate_transformation.functions.transform \
    import get_cartesian_distance
from benchmark.functions import get_nc_curvature, plot_curvature

# TODO fix get_cart_distance, does not work when import
# gives x and y in opposite shapes than what should and 
# actual distances are like 5000-9000km instead of 200km
# Transform lat, lon to be centered around the N Pole
x_N, y_N = get_cartesian_distance(lon_dom, lat_dom)

# Save .nc datasets
get_nc_curvature('spherical_surface', surface_sphere, x_N, y_N)
get_nc_curvature('spherical_ocean', ocean_sphere, x_N, y_N)
get_nc_curvature('spherical_Moho', Moho_sphere, x_N, y_N)
get_nc_curvature('spherical_bottom', bottom_sphere, x_N, y_N)

# test by plotting 
plot_curvature(lat_dom, lon_dom, bottom_sphere)
# %% Get curvature for ellipsoid
from coordinate_transformation.functions.get_spherical \
    import radius_cnt, wgs84, geographic_to_geocentric
from benchmark.functions import \
    get_curvature_wgs84

# BUG in getting cartesian distances for ellipsoid case
# BUG in curv_wgs84, don't get 0 at the centre

surface_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6370287.273) - 32.4975
ocean_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6365567.273) - 32.4734
Moho_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6358117.273) - 32.4353
bottom_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6270287.273) - 31.9873

# %% Save curvature files for oblate Earth
from netCDF4 import Dataset
import numpy as np 
import datetime as dt
from coordinate_transformation.functions.get_spherical \
    import geographic_to_geocentric, wgs84
from coordinate_transformation.functions.get_rotation \
    import get_cartesian_distance
from benchmark.functions import \
    get_nc_curvature

# Transform lat, lon to be centered around the N Pole
# x_N, y_N = get_cartesian_distance(lon_dom, lat_dom)

# Save .nc datasets
get_nc_curvature('ellipsoid_surface', surface_ellipsoid, x_N, y_N)
get_nc_curvature('ellipsoid_ocean', ocean_ellipsoid, x_N, y_N)
get_nc_curvature('ellipsoid_Moho', Moho_ellipsoid, x_N, y_N)
get_nc_curvature('ellipsoid_bottom', bottom_ellipsoid, x_N, y_N)

# %% Check netcdf variables 
from pathlib import Path
import netCDF4 as nc4 
from benchmark.dom_calculations.functions import \
    check_nc

check_nc('benchmark/input_files/relabelling/', \
    'spherical_surface.nc')
check_nc('benchmark/input_files/relabelling/', \
    'ellipsoid_surface.nc')

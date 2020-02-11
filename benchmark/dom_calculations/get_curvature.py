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
lat_max = 39.5
lat_min = 35.5
lon_max = -13.8
lon_min = -19.2
bounds = [lat_max, lat_min, lon_max, lon_min]

(lat_dom, lon_dom, bathy_dom) = truncate_domain(lat_Prt, \
    lon_Prt, bathy_Prt, bounds)

# %% Get arrays with distances to the curve for sphere
import numpy as np 
from math import pi
from benchmark.dom_calculations.functions import \
    get_curvature

surface_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6370.107295)
ocean_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6365.387295)
Moho_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6357.937295)
bottom_sphere = get_curvature(lat_dom, \
    lon_dom, radius = 6270.107295)

# %% Save curvature values for the spherical case
from netCDF4 import Dataset
import numpy as np 
import datetime as dt
from coordinate_transformation.functions.get_spherical \
    import geographic_to_geocentric, wgs84
from coordinate_transformation.functions.get_rotation \
    import get_cartesian_distance
from benchmark.dom_calculations.functions import \
    get_nc_curvature

# Transform lat, lon to be centered around the N Pole
x_N, y_N = get_cartesian_distance(lon_dom, lat_dom)

# Save .nc datasets
get_nc_curvature('spherical_surface', surface_sphere)
get_nc_curvature('spherical_ocean', ocean_sphere)
get_nc_curvature('spherical_Moho', Moho_sphere)
get_nc_curvature('spherical_bottom', bottom_sphere)

# %% Get curvature for ellipsoid
from coordinate_transformation.functions.get_spherical \
    import radius_cnt, wgs84, geographic_to_geocentric
from benchmark.dom_calculations.functions import \
    get_curvature_wgs84

# BUG in getting cartesian distances for ellipsoid case

surface_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6370.287273)
ocean_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6365.387295)
Moho_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6357.937295)
bottom_ellipsoid = get_curvature_wgs84(lat_dom, \
    lon_dom, radius = 6270.107295)

# %% Save curvature files for oblate Earth
from netCDF4 import Dataset
import numpy as np 
import datetime as dt
from coordinate_transformation.functions.get_spherical \
    import geographic_to_geocentric, wgs84
from coordinate_transformation.functions.get_rotation \
    import get_cartesian_distance
from benchmark.dom_calculations.functions import \
    get_nc_curvature

# Transform lat, lon to be centered around the N Pole
# x_N, y_N = get_cartesian_distance(lon_dom, lat_dom)

# Save .nc datasets
get_nc_curvature('ellipsoid_surface', surface_ellipsoid)
get_nc_curvature('ellipsoid_ocean', ocean_ellipsoid)
get_nc_curvature('ellipsoid_Moho', Moho_ellipsoid)
get_nc_curvature('ellipsoid_bottom', bottom_ellipsoid)

# %% Check netcdf variables 
from pathlib import Path
import netCDF4 as nc4 
from benchmark.dom_calculations.functions import \
    check_nc

check_nc('benchmark/input_files/relabelling/', \
    'spherical_surface.nc')
check_nc('benchmark/input_files/relabelling/', \
    'ellipsoid_surface.nc')

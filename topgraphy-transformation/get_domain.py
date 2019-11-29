from pathlib import Path
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma

#%%  Import global elevation files using Dataset
data_folder = Path('/Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/top_data/GEBCO_2019/')
file2open = data_folder / 'GEBCO_2019.nc' #file with location
nc_GEBCO = nc4.Dataset(file2open, 'r')
raw_lat = nc_GEBCO.variables['lat'][:] # in degrees N
raw_lon = nc_GEBCO.variables['lon'] [:]# in degrees E 
raw_elevation = nc_GEBCO.variables['elevation'] [:] 
# in m as height above reference ellipsoid

# Unmask arrays
raw_lat = np.ma.getdata(raw_lat)
raw_lon = np.ma.getdata(raw_lon)
raw_elevation = np.ma.getdata(raw_elevation)

#%% Truncate domain
# Define domain
lat_max = 41.4 
lat_min = 35.5
lon_max = -14.5
lon_min = -22
bounds = [lat_max, lat_min, lon_max, lon_min]

def find_nearest(array, value):
    #Find index of the element nearest to the value specified
    array = np.asarray(array)
    indx = (np.abs(array - value)).argmin()
    return indx

def truncate_domain(lat, lon, value, bounds):
    """ 
    Input lat, lon, elevatin value of interest (topography, depth to Moho 
    etc.), and domain boundaries as an array of 4 elements, in the order: 
    max latitude, min latitude, max longitude, min longitude. 
    """
    #Find indices
    indx_lat_max = find_nearest(lat, bounds[0])
    indx_lat_min = find_nearest(lat, bounds[1])
    indx_lon_max = find_nearest(lon, bounds[2])
    indx_lon_min = find_nearest(lon, bounds[3])
    #Truncate domain
    lat_domain = lat[indx_lat_min:indx_lat_max]
    lon_domain = lon[indx_lon_min:indx_lon_max]
    n = len(value)
    if n == len(lat):
        value_domain = value[indx_lat_min:indx_lat_max, indx_lon_min:indx_lon_max]
    elif n == len(lon): 
        value_domain = value[indx_lon_min:indx_lon_max, indx_lat_min:indx_lat_max]
    else: 
        print('Array must have same dimensions as lat and lon')
    return(lat_domain, lon_domain, value_domain)

#Truncate domain
(lat_Prt, lon_Prt, bathy_Prt) = truncate_domain(raw_lat, raw_lon, raw_elevation, bounds)

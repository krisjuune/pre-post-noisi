from pathlib import Path
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma

def find_nearest(array, value):
    """
    Find index of the element in a one-dimensional array
    nearest to the value specified. 
    """
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

def relative_depth(elevation, reference_value): 
    """ 
    Return depths relative to reference value, given the 
    elevation. Elevation dataset must be negative below
    surface. Output dataset is positive downwards. 
    """
    # write code
    return(relative_depth)
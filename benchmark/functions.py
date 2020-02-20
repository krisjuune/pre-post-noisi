# %% Calculations
import numpy as np 
import numpy.ma as ma
from math import pi
from coordinate_transformation.functions.get_spherical \
    import radius_cnt, wgs84, geographic_to_geocentric
from coordinate_transformation.functions.get_domain \
    import find_nearest

def get_curvature(lat, lon, radius = 6370287.272978241, \
    theta = 37.5, phi = -16.5):
    """
    Function to calculate the curvature relative to a 
    flat surface at the given radius assuming a sphere 
    with the given radius. Inputs include arrays of 
    latitude, longitude, and a radius. Function returns 
    array of depths relative to this flat surface with 
    the dimensions of lon, lat. 
    Units in degrees for angles, distances same as radius. 
    Default radius calculated at default geographic theta.   
    """
    # preallocate output array
    curvature = np.zeros((len(lon), \
        len(lat)), float) 
    
    # transform to geocentric
    lat = geographic_to_geocentric(lat)
    theta = geographic_to_geocentric(theta)

    # convert to radians
    lon = pi/180*lon
    lat = pi/180*lat
    phi = pi/180*phi
    theta = pi/180*theta

    # loop over the lats and lons
    for i in range(len(lon)):
        for j in range(len(lat)):
            # find angle between point i,j and 
            # centre
            a = radius*np.sin(lon[i] - phi)
            b = radius*np.sin(lat[j] - theta)
            c = np.sqrt(np.square(a) + \
                np.square(b))
            # arcsin(x), x has to be [-1,1]
            alpha = np.arcsin(c/radius)
            # calculate depth to curve from flat 
            # surface
            y = radius/np.cos(alpha) - radius
            x = y*np.cos(alpha)
            curvature [i,j] = x*(-1)
    
    return(curvature)

def get_curvature_wgs84(lat, lon, radius = 6370287.272978241, \
    theta = 37.5, phi = -16.5):
    """
    Function to calculate the curvature relative to a 
    flat surface at the given radius for an ellipsoid 
    defined by wgs84. Inputs include arrays of latitude, 
    longitude, and a radius. Function returns array of 
    depths relative to this flat surface with the 
    dimensions of lon, lat. 
    Units in degrees for angles, distances same as radius. 
    Default radius calculated at default geographic theta. 
    """
    # preallocate output array
    curvature = np.zeros((len(lon), \
        len(lat)), float) 

    # transform to geocentric
    lat = geographic_to_geocentric(lat)
    theta = geographic_to_geocentric(theta)

    # convert to radians
    lon = pi/180*lon
    lat = pi/180*lat
    phi = pi/180*phi
    theta = pi/180*theta 

    # loop over the lats and lons
    for i in range(len(lon)):
        for j in range(len(lat)):
            # find radius at j-th latitude
            if round(radius/1000, 3) == radius_cnt(theta)/1000: 
                # when look at the surface curvature
                # centred around default lat, lon
                a = wgs84()[0]
                b = wgs84()[1]
            else:
                # for when looking at shallower levels
                # centred about any lat, lon
                r_theta = radius_cnt(theta)/1000
                a = wgs84()[0]*radius/r_theta
                b = wgs84()[1]*radius/r_theta

            radius_j = np.sqrt((a**2*(np.cos(lat[j])**2)) + \
                (b**2*(np.sin(lat[j])**2)))/1000

            # find angle between point i,j and centre
            l1 = abs(radius*np.tan(lon[i] - phi))
            l2 = abs(radius*np.tan(lat[j] - theta))
            l3 = np.sqrt(np.square(l1) + \
                np.square(l2))
            alpha = np.arctan(l3/radius)
            # Checked and up to alpha everything seems to 
            # be working

            # calculate depth to curve from flat surface
            y = radius/np.cos(alpha) - radius_j
            x = y*np.cos(alpha)
            curvature [i,j] = x*(-1) 
    
    # Cannot seem to find reason why curvature !=0 at 
    # (theta, phi), so just substituting that value 
    # from all elements, tested this against spherical
    # case and get same values to ~10m accuracy
    m = find_nearest(lon, phi)
    n = find_nearest(lat, theta)
    curvature = curvature - curvature [m,n]
    
    return(curvature)

# %% Save as netCDF files, add a check function
import netCDF4 as nc4 
import datetime as dt 
from mpl_toolkits.basemap import Basemap
from pathlib import Path

def get_nc_curvature(filename, curvature_variable, x_var, y_var):
    """
    Writes a netCDF4 file with x_distance, y_distance, 
    and curvature. filename should be a string (with
    no file extension) and curvature_variable an array 
    containing the calculated curvature values. 
    """
    
    # Create .nc file
    import netCDF4 as nc4
    f = nc4.Dataset(filename + '.nc','w', format = 'NETCDF4')
    f.description = 'Curvature calculated relative to tangent' + \
        ' surface at the centre of domain assuming spherical Earth'

    # Create dimensions
    f.createDimension('x', len(x_var))
    f.createDimension('y', len(y_var))

    # Create variables, 'f4' for single precision floats, i.e. 32bit
    curvature = f.createVariable('curvature', 'f4', ('x', 'y'))
    curvature [:] = curvature_variable
    x = f.createVariable('x', 'f4', 'x')
    x [:] = x_var
    y = f.createVariable('y', 'f4', 'y')
    y [:] = y_var

    # Add attributes to the file
    today = dt.datetime.now()
    f.history = "Created " + today.strftime("%d/%m/%y")
    #Add local attributes to variable instances
    curvature.units = 'm'
    x.units = 'm'
    y.units = 'm'

    f.close()

def check_nc(path, filename):
    from pathlib import Path
    path = Path(path)
    f = nc4.Dataset(path / filename, 'r')
    for i in f.variables:
        print(i, f.variables[i].units, \
            f.variables[i].shape)

# %% Plotting 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def plot_geographic(lat, lon, data, filename, \
    lat_max = 39.5, lat_min = 35.5, lon_max = -14, \
    lon_min = -19, cbar_label = 'Bathymetry (km)'):
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
    cmap = 'viridis'
    Lon, Lat = np.meshgrid(lon, lat)
    fig.pcolormesh(Lon, Lat, data, latlon = True, \
        cmap = cmap)
    # Colorbar construction
    i = ax.imshow(data, interpolation='nearest')
    cbar = fig.colorbar(i, shrink = 0.5, aspect = 5)
    cbar.set_label(cbar_label, rotation = 270, labelpad=15, y=0.45)

    plt.savefig(filename, dpi = 600)
    plt.show()

# TODO fix these dependencies issues, had to run to define each of 
# the functions below in order to be able to use the plotting function

def wgs84(): 
    """
    WGSS84 coordinate system with Greenwich as lon = 0.
    Define Earth's semi-major, semi-minor axes, and
    inverse flattening, eccentricity in this order. 
    """
    # set semi-major axis of the oblate spheroid Earth, in m
    a = 6378137.0
    # set semi-minor axis of the oblate spheroid Earth, in m
    b = 6356752.314245
    # calculate inverse flattening f
    f = a/(a-b)
    # calculate squared eccentricity e
    e_2 = (a**2-b**2)/a**2
    return(a,b,e_2,f)

def geographic_to_geocentric(lat):
    """
    Calculate latitude defined in the wgs84 coordinate 
    system given the geographic latitude. Input and 
    output latitude in degrees. 
    """
    e_2 = wgs84()[2] # eccentricity as defined by wgs84 
    lat = np.rad2deg(np.arctan((1 - e_2) * np.tan(np.deg2rad(lat))))
    return lat

def radius_cnt(lat):
    """
    Get radius at latitude lat for the Earth as defined 
    by the wgs84 system. 
    """
    a = wgs84()[0]
    b = wgs84()[1]
    # Calculate radius for reference ellipsoid, in m 
    lat = pi/180*lat
    # for i in range(len(lat)): 
    r_cnt = np.sqrt((a**2*(np.cos(lat)**2)) + \
        (b**2*(np.sin(lat)**2)))
    return(r_cnt)

def get_cartesian_distance(lon, lat, \
    src_lat = 37.5, src_lon = -16.5):
    """
    Calculate distance of each point of lat and lon
    from the source location on a flat surface, 
    tangential to the source. Returns x (lon), y 
    (lat) in km for AxiSEMCartesian. 
    """
    # transform to geocentric
    lat = geographic_to_geocentric(lat)
    src_lat = geographic_to_geocentric(src_lat)

    # find radius at source
    r_greatcircle = radius_cnt(src_lat)/1000
    # find radius of small circle at source lat
    r_smallcircle = r_greatcircle*np.cos(np.deg2rad(src_lat))
    # convert differences in angles to radians
    phi = pi/180*lon - pi/180*src_lon
    theta = pi/180*lat - pi/180*src_lat

    # preallocate output arrays
    x = np.zeros(len(phi), float)
    y = np.zeros(len(theta), float)
    # find distances
    x = r_smallcircle*np.tan(phi)
    y = r_greatcircle*np.tan(theta)

    return(x,y)

def plot_curvature(lat, lon, curvature, src_lat = 37.5, \
    src_lon = -16.5, cbar_label = 'Curvature (m)', \
    filename = 'noname'):
    """
    Function to plot a 3d surface once transformed 
    into Cartesian distances. Figure saved as png if 
    filename is not noname. BUG fixed with transposing. 
    """
    # TODO error with get_cart_dist so just added it in here
    # Transform lat, lon to be centered around the N Pole
    (x, y) = get_cartesian_distance(lon, lat, \
        src_lat, src_lon)
    x, y = np.meshgrid(x, y)
    x = x.transpose()
    y = y.transpose()

    # Create figure handle
    fig = plt.figure()
    ax = plt.gca(projection = '3d')
    # TODO how to scale the axes, so z not so exaggerated
    # Plot 
    surf = ax.plot_surface(x, y, curvature, \
        cmap = 'viridis')
    # Add colorbar
    cbar = fig.colorbar(surf, shrink = 0.5, aspect = 5)
    cbar.set_label(cbar_label, rotation = 270, labelpad=15, \
        y=0.45)

    plt.show()

    if filename != 'noname':
        plt.savefig((filename + '.png'), dpi = 600)
    
# %% Processing functions
import numpy as np 
from pathlib import Path
from math import pi
from coordinate_transformation.functions.get_spherical \
    import wgs84

def station_data(path, station):
    """
    Function that retrieves the seismic data from station
    'station' given the relative path 'path', both inputs 
    are strings. This works for II type (not IU) stations.
    Function returns the data array. 
    """
    path = Path(path)
    file = 'II.' + station + '.RTZ.ascii'
    # file handle
    file = path/file 
    # Open file and retrieve data
    raw_data = open(file, 'r')
    raw_data = raw_data.read()
    raw_data = np.fromstring(raw_data, dtype = float, sep=' ')
    # nr of columns is always 4 since time, rr, tt, zz
    m = int(len(raw_data)/4) 
    # preallocate output array
    data = np.zeros(((m),4), float)
    # retrieve data (which has been sorted row-wise)
    # and sort it column-wise, returning every 4th element
    data[:,0] = raw_data[0::4] 
    data[:,1] = raw_data[1::4]
    data[:,2] = raw_data[2::4]
    data[:,3] = raw_data[3::4]
    return(data)

# Calculate the length of one degree of lat and lon as a function of lat
def len_deg_lon(lat):
    """
    Calculates length of one degree of longitude
    at latitudes lat. Input lat must be an array
    of integers. 
    """
    e_2 = wgs84()[2]
    a = wgs84() [0]
    # This is the length of one degree of longitude 
    # approx. after WGS84, at latitude lat
    # in m
    lat = pi/180*lat
    dlon = (pi*a*np.cos(lat))/180*np.sqrt((1-e_2*np.sin(lat)**2))
    return np.round(dlon,5)

def len_deg_lat(lat):
    """
    Calculates length of one degree of latitude
    at latitudes lat. Input lat must be an array
    of integers. 
    """
    # This is the length of one degree of latitude 
    # approx. after WGS84, between lat-0.5deg and lat+0.5 deg
    # in m
    lat = pi/180*lat
    dlat = 111132.954 - 559.822 * np.cos(2*lat) + 1.175*np.cos(4*lat)
    return np.round(dlat,5)

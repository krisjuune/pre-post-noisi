# %% Calculations
import numpy as np 
import numpy.ma as ma
from math import pi
from coordinate_transformation.functions.get_spherical \
    import radius_cnt, wgs84, geographic_to_geocentric
from coordinate_transformation.functions.get_domain \
    import find_nearest

def get_curvature(lat, lon, radius = 6370.287272978241, \
    theta = 37.5, phi = -16.5):
    """
    Function to calculate the curvature 
    relative to a flat surface at the given 
    radius assuming a sphere with the given 
    radius. Inputs include arrays of latitude, 
    longitude, and radius. Function returns 
    array of depths relative to this flat 
    surface with the dimensions of lon, lat.
    Units in degrees for angles, km for 
    distances.  
    """
    # preallocate output array
    curvature = np.zeros((len(lon), \
        len(lat)), float) 
    
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
            curvature [i,j] = x 
    
    return(curvature)

def get_curvature_wgs84(lat, lon, radius = 6370.287272978241, \
    theta = 37.5, phi = -16.5):
    """
    Function to calculate the curvature relative to 
    a flat surface at the given radius for an 
    ellipsoid defined by wgs84. Inputs include 
    arrays of latitude, longitude, and a radius. 
    Function returns array of depths relative to 
    this flat surface with the dimensions of lon, 
    lat. 
    Units in degrees for angles, km for distances. 
    Default radius for latitude theta and longitude 
    phi. 
    """
    # preallocate output array
    curvature = np.zeros((len(lon), \
        len(lat)), float) 

    # transform to geocentric coordinates
    lat = geographic_to_geocentric(lat)

    # convert to radians
    lon = pi/180*lon
    lat = pi/180*lat
    phi = pi/180*phi
    theta = pi/180*theta 

    # loop over the lats and lons
    for i in range(len(lon)):
        for j in range(len(lat)):
            # find radius at j-th latitude
            if round(radius, 3) == radius_cnt(theta)/1000: 
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
            curvature [i,j] = x 
    
    # Cannot seem to find reason why curvature !=0 at 
    # (theta, phi), so just substituting that value 
    # from all elements, tested this against spherical
    # case and get same values to ~10m accuracy
    m = find_nearest(lon, phi)
    n = find_nearest(lat, theta)
    curvature = curvature - curvature [m,n]
    
    return(curvature)

# %% Save as netCDF files, add a check function

def get_nc_curvature(filename, curvature_variable):
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
    f.createDimension('x', len(x_N))
    f.createDimension('y', len(y_N))

    # Create variables, 'f4' for single precision floats, i.e. 32bit
    curvature = f.createVariable('curvature', 'f4', ('x', 'y'))
    curvature [:] = curvature_variable
    x = f.createVariable('x', 'f4', 'x')
    x [:] = x_N
    y = f.createVariable('y', 'f4', 'y')
    y [:] = y_N

    # Add attributes to the file
    today = dt.datetime.now()
    f.history = "Created " + today.strftime("%d/%m/%y")
    #Add local attributes to variable instances
    curvature.units = 'km'
    x.units = 'km'
    y.units = 'km'

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

def plot_curvature(lat, lon, curvature, src_lat = 37.5, \
    src_lon = -16.5, cbar_label = 'Curvature (km)', \
    filename = 'noname'):
    """
    Function to plot a 3d surface once transformed 
    into Cartesian distances. Figure saved as png 
    if filename is not noname.  
    """
    
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
    surf = ax.plot_surface(x, y, curvature*(-1), \
        cmap = 'viridis')
    # Add colorbar
    cbar = fig.colorbar(surf, shrink = 0.5, aspect = 5)
    cbar.set_label(cbar_label, rotation = 270, labelpad=15, \
        y=0.45)

    plt.show()

    if filename != 'noname':
        plt.savefig((filename + '.png'), dpi = 600)
    
    

# Import lat and lon files of the domain off Portugal

import numpy as np
from pathlib import Path
from scipy.io import netcdf
# from netCDF4 import Dataset

# Read GEBCO_2019.nc netCDF4 using Dataset 
data_folder = Path('/Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/top_data/GEBCO_2019/')
file2open = data_folder / 'GEBCO_2019.nc' #file with location
nc_GEBCO = Dataset(file2open, 'r')
raw_lat = nc_GEBCO.variables['lat'][:] # import lat 
raw_lon = nc_GEBCO.variables['lon'] # import lat 
raw_elevation = nc_GEBCO.variables['elevation'] # import lat 

def var_read():
    f = open('lat35.5_41.4N.txt', 'r')
    # if f.mode == "r" #check if in default read mode
    # use read() to read the contents of the file
    lat_Prt = f.read()
    check = len(lat_Prt)
    print(check)

lat_test = [line.rstrip('\n') for line in open('lat_Prt.txt')]

# def var_read():
#     f = open('lon14.5_22W.txt', 'r')
#     # if f.mode == "r" #check if in default read mode
#     # use read() to read the contents of the file
#     lon_Prt = f.read()
#     check = len(lon_Prt)
#     print(check)

# if __name__ = '__main__'
# main()

# print("before import")
# import math

# print("before functionA")
# def functionA():
#     print("Function A")

# print("before functionB")
# def functionB():
#     print("Function B {}".format(math.sqrt(100)))

# print("before __name__ guard")
# if __name__ == '__main__':
#     functionA()
#     functionB()
# print("after __name__ guard")





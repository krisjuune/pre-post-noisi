from pathlib import Path
from netCDF4 import Dataset

# Read GEBCO_2019.nc netCDF4 using Dataset 
data_folder = Path('/Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/top_data/GEBCO_2019/')
file2open = data_folder / 'GEBCO_2019.nc' #file with location
nc_GEBCO = Dataset(file2open, 'r')
raw_lat = nc_GEBCO.variables['lat'][:] # import lat 
raw_lon = nc_GEBCO.variables['lon'] # import lat 
raw_elevation = nc_GEBCO.variables['elevation'] # import lat 

raw_elevation.shape #Check array dimensions to see if reading .nc worked, 43200 (lat) x 86400 (lon)

lat_Prt = raw_lat[30120:31536] #Indices for 35.5-41.4N obtained in MatLAB
lon_Prt = raw_lon[37920:39720] #Indices for -22 - -14.5E obtained in MatLAB
bathy_Prt = raw_elevation[30120:31536, 37920:39720]

bathy_Prt.shape #check array dimensions, 1416 x 1800


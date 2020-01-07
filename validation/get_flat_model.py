# %% Get average bathymetry and depth to Moho 
import numpy as np
import pickle
from pathlib import Path
# Could weight the average calculations by the len of lat and lon at each point 
# But do it here in a rubbish manner by just taking the unweighted average

# Get average bathymetry for validation
with open('coordinate_transformation/bathy_Prt', 'rb') as f:
    bathy_Prt = pickle.load(f)
bathy_Prt = np.ma.getdata(bathy_Prt)
print(np.mean(bathy_Prt)*(-1), 'is the average depth in m')

# Get average Moho depth for validation
data_folder = Path('coordinate_transformation/raw_data/crust1.0/')
moho_file = data_folder / 'depthtomoho.xyz'
raw_moho = open(moho_file, 'r')
raw_moho = np.ma.getdata(raw_moho) 
# in columns lat, lon, depth

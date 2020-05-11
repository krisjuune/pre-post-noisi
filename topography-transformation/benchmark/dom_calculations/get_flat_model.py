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

# %% Get 3layers model properties

def get_3layers(layer1, layer2, Earth_radius):
    """ 
    Function that prints the properties to enter into the .bm
    model file, given the depths of the layers from the Earth's
    surface and the Earth radius. Inputs must be floats. 
    """
    if type(layer1) == float \
        and type(layer2) == float \
        and type(Earth_radius) == float:
        l0 = str(Earth_radius) + ' '
        l1 = str(Earth_radius - layer1) + ' '
        l2 = str(Earth_radius - layer2) + ' '
        l3 = str(0.) + ' '
        print('Discontinuities are at: ' \
            + l0 + l1 + l2 + l3)
    else:
        print('Error, inputs must be floats')

r_Earth = 6371.0 # km
# with open('coordinate_transformation/r38', 'rb') as f:
#     r38 = pickle.load(f)
# r38 = np.ma.getdata(r38)
av_bathy = 4.72 # km, average ocean depth in dom
av_Moho = 12.17 # km, average Moho depth in dom

get_3layers(av_bathy, av_Moho, r_Earth)

# %%

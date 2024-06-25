"""
Author : José CUNHA TEIXEIRA
License : SNCF Réseau, UMR 7619 METIS
Date : November 30, 2023
"""

import os
import numpy as np
from obspy import read
from folders import data_dir, results_dir
from display import display_seismic_wiggle_fromStream



### ACTIVE ----------
# profile = "P1"
# record = "1001"

# stream = read(f"{data_dir}/{profile}/Active/{record}.dat")
# offsets = np.arange(0.0, 24, 0.25)
# path = f"{results_dir}/{profile}/Active/{record}/"
# if not os.path.exists(path):
#     os.makedirs(path)
# display_seismic_wiggle_fromStream(stream, offsets, norm_method='trace', display_method='save', scale=1, path=f"{path}/{record}_VSG.png")



### PASSIVE ----------
profile = "P1"


# To display a VSG
params = "all-Sl5-Ss1-FK30-PWS2"
N_MASW = 96
path = f"{results_dir}/{profile}/Passive/{params}/W{N_MASW}/"
name = f"{profile}_VSG"
stream = read(f"{results_dir}/{profile}/Passive/{params}/W{N_MASW}/{name}.segy")

## To display a RAW Record
# record = "101"
# stream = read(f"{data_dir}/{profile}/Passive/{record}.dat")
# path = f"{results_dir}/{profile}/Passive/{record}/"
# if not os.path.exists(path):
#     os.makedirs(path)
# name = f"{record}_SG"


offsets = np.arange(0.0, 24, 0.25)
display_seismic_wiggle_fromStream(stream, offsets, norm_method='trace', display_method='save', scale=1, path=f"{path}/{name}.png")
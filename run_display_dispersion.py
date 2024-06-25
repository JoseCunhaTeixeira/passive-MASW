"""
Author : José CUNHA TEIXEIRA
License : SNCF Réseau, UMR 7619 METIS
Date : November 30, 2023
"""

import os
import numpy as np
from folders import results_dir
from display import display_dispersion_img



### ACTIVE ----------
# profile = "P1"
# N_MASW = 96

# dir = f"{results_dir}/P1/Active/W{N_MASW}/"
# FV = np.load(dir + "P1_dispIm.npy")
# fs = np.load(dir+ "P1_fs.npy")
# vs = np.load(dir+ "P1_vs.npy")
# mode_files = [file for file in os.listdir(dir) if file.startswith(f'{profile}_M') and not file.endswith('wlgt.pvc')]
# modes = []
# for mode in mode_files:
#     modes.append(np.loadtxt(dir + mode))
# display_dispersion_img(FV, fs, vs, *modes, display_method='save', path=dir+"P1_dispIm_picked.png", normalization="Frequency", errbars=False)
# display_dispersion_img(FV, fs, vs, *modes, display_method='save', path=dir+"P1_dispIm_picked_err.png", normalization="Frequency", errbars=True)
# display_dispersion_img(FV, fs, vs, display_method='save', path=dir+"P1_dispIm.png", normalization="Frequency", errbars=False)



### PASSIVE ----------
profile = "P1"
params = "all-Sl5-Ss1-FK30-PWS2"
N_MASW = 96

dir = f"{results_dir}/{profile}/Passive/{params}/W{N_MASW}/"
FV = np.load(dir + "P1_dispIm.npy")
fs = np.load(dir + "P1_fs.npy")
vs = np.load(dir + "P1_vs.npy")
mode_files = [file for file in os.listdir(dir) if file.startswith(f'{profile}_M') and not file.endswith('wlgt.pvc')]
modes = []
for mode in mode_files:
    modes.append(np.loadtxt(dir + mode))
display_dispersion_img(FV, fs, vs, *modes, display_method='save', path=dir+"P1_dispIm_picked.png", normalization="Frequency", errbars=False)
display_dispersion_img(FV, fs, vs, *modes, display_method='save', path=dir+"P1_dispIm_picked_err.png", normalization="Frequency", errbars=True)
display_dispersion_img(FV, fs, vs, display_method='save', path=dir+"P1_dispIm.png", normalization="Frequency", errbars=False)
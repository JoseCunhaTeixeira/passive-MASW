"""
Author : José CUNHA TEIXEIRA
License : SNCF Réseau, UMR 7619 METIS
Date : November 30, 2023
"""

import matplotlib.pyplot as plt
import numpy as np
from display import display_dispersion_img

plt.rcParams.update({'font.size': 22})
colors = ["g", "r", "b", "m", "y", "w", "c", "o"]


### ACTIVE ----------
dir = "/.../data_passive-MASW/Results/P1/Active/W96/"
FV = np.load(dir + "P1_dispIm.npy")
fs = np.load(dir+ "P1_fs.npy")
vs = np.load(dir+ "P1_vs.npy")
mode0 = np.loadtxt(dir + 'P1_M0.pvc')
mode1 = np.loadtxt(dir + 'P1_M1.pvc')
mode2 = np.loadtxt(dir + 'P1_M2.pvc')
display_dispersion_img(FV, fs, vs, mode0, mode1, display_method='save', path=dir+"P1_dispIm_picked.png", normalization="Frequency", errbars=False)
display_dispersion_img(FV, fs, vs, mode0, mode1, display_method='save', path=dir+"P1_dispIm_picked_err.png", normalization="Frequency", errbars=True)
display_dispersion_img(FV, fs, vs, display_method='save', path=dir+"P1_dispIm.png", normalization="Frequency", errbars=False)



### PASSIVE ----------
# dir = "/.../data_passive-MASW/Results/P1/Passive/all-Sl5-Ss1-FK30-PWS2/W96/"
# FV = np.load(dir + "P1_dispIm.npy")
# fs = np.load(dir + "P1_fs.npy")
# vs = np.load(dir + "P1_vs.npy")
# mode0 = np.loadtxt(dir + 'P1_M0.pvc')
# mode1 = np.loadtxt(dir + 'P1_M1.pvc')
# mode2 = np.loadtxt(dir + 'P1_M2.pvc')
# mode3 = np.loadtxt(dir + 'P1_M3.pvc')

# display_dispersion_img(FV, fs, vs, mode0, display_method='save', path=dir+"P1_dispIm_picked.png", normalization="Frequency", errbars=False)
# display_dispersion_img(FV, fs, vs, mode0, display_method='save', path=dir+"P1_dispIm_picked_err.png", normalization="Frequency", errbars=True)
# display_dispersion_img(FV, fs, vs, display_method='save', path=dir+"P1_dispIm.png", normalization="Frequency", errbars=False)



### ADDITION ACTIVE DIRECT AND INVERSE SHOTS ----------
# a_dir = "/.../data_passive-MASW/Results/P1/Active/1001/W96/"
# a_dir_FV = np.load(a_dir + "1001_dispIm.npy")
# a_dir_fs = np.load(a_dir+ "1001_fs.npy")
# a_dir_vs = np.load(a_dir+ "1001_vs.npy")

# a_inv = "/.../data_passive-MASW/Results/P1/Active/1002/W96/"
# a_inv_FV = np.load(a_inv + "1002_dispIm.npy")
# a_inv_fs = np.load(a_inv+ "1002_fs.npy")
# a_inv_vs = np.load(a_inv+ "1002_vs.npy")

# a_FV = (a_dir_FV + a_inv_FV)/2

# path = "/.../data_passive-MASW/Results/P1/Active/W96/"
# np.save(path+"P1_dispIm.npy", a_FV)
# np.save(path+"P1_fs.npy", a_dir_fs)
# np.save(path+"P1_vs.npy", a_dir_vs)
# display_dispersion_img(a_FV, a_dir_fs, a_dir_vs, display_method='save', normalization="Frequency", path=path+"P1_dispIm.png")
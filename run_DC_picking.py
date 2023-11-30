import csv
import sys

import colormaps as cmaps
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('./../Tools/')
from signal_processing import pick, resamp

plt.rcParams.update({'font.size': 28})


dir = "/.../project_name/Results/1/Passive/all-Sl5-Ss1-FK30-Nsegsmax50-PWS2/W96_0-95/"
profile = "P1"
mode = 0 #mode number

FV = np.load(dir + f"{profile}_dispIm.npy")
fs = np.load(dir+ f"{profile}_fs.npy")
vs = np.load(dir+ f"{profile}_vs.npy")


dx = 0.25 #space between geophones
Nx = 96 #number of geophones



### One curve picking -----------
for i, col in enumerate(FV.T):
    FV[:, i] = FV[:, i] / np.max(col)
    
f_picked, v_picked, err = pick(FV, fs, vs, dx, Nx, smooth=False)

with open(dir + f'{profile}_M{mode}.pvc', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(f_picked, v_picked, err))

f_resamp, v_resamp, err_resamp = resamp(f_picked, v_picked, err, type='wavelength-log')
with open(dir + f'{profile}_M{mode}_log-wlgt.pvc', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(f_resamp, v_resamp, err_resamp))

f_resamp, v_resamp, err_resamp = resamp(f_picked, v_picked, err, type='wavelength')
with open(dir + f'{profile}_M{mode}_wlgt.pvc', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(f_resamp, v_resamp, err_resamp))


plt.figure(figsize=(16, 9), dpi=100)
extent = [fs[0], fs[-1], vs[0], vs[-1]]
plt.imshow(np.flipud(FV**2), cmap=cmaps.wh_bl_gr_ye_re, extent=extent, aspect="auto")
# plt.plot(f_picked, v_picked, "+w")
plt.errorbar(f_picked, v_picked, err, fmt="+w", ecolor='black', elinewidth=0.5, mew=0.5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase velocity (m/s)")
plt.title(f"Dispersion image with picked curve")
plt.colorbar()
plt.xticks(np.arange(min(fs), max(fs)+20.0, 20.0))
v_ticks = np.arange(0, max(vs), 200.0)
v_ticks[0] = 1
plt.yticks(v_ticks)
plt.tight_layout()
plt.show()



### MASW profile -----------
# for i in np.arange(0, 91, 1):
#     print(f'{int(i)}-{int(i)+5}')
#     dir = f"/home/jteixeira/Documents/PhD_Monitoring/Work/Processing/Data/Sarrebourg/2023-04-03_Sarrebourg/Results/P1/Passive/all-Sl5-Ss1-FK0-Nsegsmax20-PWS2/W6_{int(i)}-{int(i)+5}/"

#     profile = "P1"
#     mode = 0

#     FV = np.load(dir + f"{profile}_dispIm.npy")
#     fs = np.load(dir+ f"{profile}_fs.npy")
#     vs = np.load(dir+ f"{profile}_vs.npy")


#     dx = 0.25
#     Nx = 96


#     for i, col in enumerate(FV.T):
#         FV[:, i] = FV[:, i] / np.max(col)
        
#     f_picked, v_picked, err = pick(FV, fs, vs, dx, Nx, smooth=False)

#     with open(dir + f'{profile}_M{mode}.pvc', 'w') as f:
#         writer = csv.writer(f, delimiter='\t')
#         writer.writerows(zip(f_picked, v_picked, err))

#     f_resamp, v_resamp, err_resamp = resamp(f_picked, v_picked, err, type='wavelength-log')
#     with open(dir + f'{profile}_M{mode}_log-wlgt.pvc', 'w') as f:
#         writer = csv.writer(f, delimiter='\t')
#         writer.writerows(zip(f_resamp, v_resamp, err_resamp))

#     f_resamp, v_resamp, err_resamp = resamp(f_picked, v_picked, err, type='wavelength')
#     with open(dir + f'{profile}_M{mode}_wlgt.pvc', 'w') as f:
#         writer = csv.writer(f, delimiter='\t')
#         writer.writerows(zip(f_resamp, v_resamp, err_resamp))


#     plt.figure(figsize=(16, 9), dpi=100)
#     extent = [fs[0], fs[-1], vs[0], vs[-1]]
#     plt.imshow(np.flipud(FV**2), cmap=cmaps.wh_bl_gr_ye_re, extent=extent, aspect="auto")
#     # plt.plot(f_picked, v_picked, "+w")
#     plt.errorbar(f_picked, v_picked, err, fmt="+w", ecolor='black', elinewidth=0.5, mew=0.5)
#     plt.xlabel("Frequency (Hz)")
#     plt.ylabel("Phase velocity (m/s)")
#     plt.title(f"Dispersion image with picked curve")
#     plt.colorbar()
#     plt.xticks(np.arange(min(fs), max(fs)+20.0, 20.0))
#     v_ticks = np.arange(0, max(vs), 200.0)
#     v_ticks[0] = 1
#     plt.yticks(v_ticks)
#     plt.tight_layout()
#     plt.show()
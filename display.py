import numpy as np
import matplotlib.pyplot as plt
from misc import verify_expected
import colormaps as cmaps

plt.rcParams.update({'font.size': 22})

CRED = "\033[91m"
CYEL = "\033[93m"
CGRE = "\033[92m"
BOLD = "\033[1m"
CEND = "\033[0m"
colors = ["g", "r", "b", "m", "y", "w", "c", "o"]
_DPI = 300




#----------------------------------------------------------------------------------------------------
def diag_print(case, str1, str2):
    if case in ("Error", "error", "ERROR"):
        return print(BOLD + CRED + "ERROR     | " + str1 + "\n          | " + str2 + "\n" + CEND)
    elif case in ("Warning", "warning", "WARNING"):
        return print(CYEL + "WARNING   | " + str1 + "\n          | " + str2 + "\n" + CEND)
    elif case in ("Info", "info", "INFO"):
        return print(CGRE + "INFO      | " + str1 + "\n          | " + str2 + "\n" + CEND)




#----------------------------------------------------------------------------------------------------
def display_seismic_wiggle_fromStream(stream, x_sensors, display_method="disp", path=None, scale=1.0, norm_method=None, **kwargs):
    for i, tr in enumerate(stream) :
        tr.stats.distance = x_sensors[i]*1000
    fig = plt.figure(figsize=(16, 9), dpi=_DPI)
    if len(kwargs) == 0:
        stream.plot(type='section', time_down=True, scale=scale, fig=fig, fillcolors=("black",None), norm_method=norm_method)
    if len(kwargs) == 2:
        verify_expected(kwargs, ["recordstart", "recordlength"])
        recordstart=kwargs["recordstart"]
        recordlength=kwargs["recordlength"]
        stream.plot(type='section', time_down=True, scale=scale, fig=fig, fillcolors=("black",None), norm_method=norm_method, recordstart=recordstart, recordlength=recordlength)
    plt.tight_layout()
    plt.xlabel("Position (m)")
    plt.ylabel("Time (s)")

    if display_method == "disp":
        plt.show()
    elif display_method == "save" and path != None:
        plt.savefig(path, format='png', dpi='figure', bbox_inches='tight')
        plt.close()




#----------------------------------------------------------------------------------------------------
def display_spectrum_img_fromArray(array, dt, x_sensors, display_method="disp", path1=None, path2=None, norm_method=None):
    Nt, Nx = array.shape
    SP = np.fft.rfft(array, axis=0)
    fs = np.fft.rfftfreq(Nt, dt)

    if norm_method == "stream":
        SP = np.abs(SP)/np.max(np.abs(SP))
    elif norm_method == "trace":
        for i, col in enumerate(SP.T):
            SP[:,i] = np.abs(col) / np.max(np.abs(col))

    extent = [x_sensors[0], x_sensors[-1], fs[0], fs[-1]]
    
    plt.figure(figsize=(16, 9), dpi=_DPI)
    plt.imshow(np.flipud(np.abs(SP)), cmap='gray_r', extent=extent, aspect="auto")
    plt.xlabel("Position (m)")
    plt.ylabel("Frequency (Hz)")
    plt.colorbar()
    plt.tight_layout()
    if display_method == "save":
        plt.savefig(path1, format='png', dpi='figure', bbox_inches='tight')
        plt.close()

    plt.figure(figsize=(16, 9), dpi=_DPI)
    plt.plot(fs, abs(SP[:,0]), fs, abs(SP[:, Nx//2]), fs, abs(SP[:,-1]), linewidth=0.5)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude")
    plt.legend(["Trace 1", f"Trace {Nx//2}", f"Trace {Nx}"])
    if display_method == "save":
        plt.savefig(path2, format='png', dpi='figure', bbox_inches='tight')
        plt.close()

    if display_method == "disp":
        plt.show()




#----------------------------------------------------------------------------------------------------
def display_dispersion_img(FV_arr, fs, vs, *args, display_method="disp", path=None, normalization=None, errbars=True):
    """
    args = mode0, mode1, ...
    """
    FV = np.copy(FV_arr)
    if normalization == "Frequency":
        for i, col in enumerate(FV.T):
            FV[:, i] = FV[:, i] / np.max(col)
    elif normalization == "Global":
        FV = FV / np.max(FV)
    plt.figure(figsize=(16, 9), dpi=_DPI)
    extent = [fs[0], fs[-1], vs[0], vs[-1]]
    plt.imshow(np.flipud(FV**2), cmap=cmaps.wh_bl_gr_ye_re, extent=extent, aspect="auto")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Phase velocity (m/s)")
    plt.xticks(np.arange(min(fs), max(fs)+50.0, 50.0))
    v_ticks = np.arange(0, max(vs), 200.0)
    v_ticks[0] = 1
    plt.yticks(v_ticks)
    plt.ylim([1, max(vs)])
    plt.tight_layout()
    if len(args) == 0 :
        if display_method == "disp":
            plt.show()
        elif display_method == "save":
            plt.savefig(path, format='png', dpi='figure', bbox_inches='tight')
            plt.close()
    else :
        for j, mode in enumerate(args):
            fs_mode = np.empty(len(mode))
            vs_mode = np.empty(len(mode))
            dc_mode = np.empty(len(mode))
            for i, line in enumerate(mode):
                fs_mode[i] = mode[i][0]
                vs_mode[i] = mode[i][1]
                dc_mode[i] = mode[i][2]
            indexes = np.where(fs_mode%np.round(fs_mode)==0)[0]
            fs_mode = fs_mode[indexes]
            vs_mode = vs_mode[indexes]
            dc_mode = dc_mode[indexes]
            if errbars == True:
                plt.errorbar(fs_mode[::1], vs_mode[::1], dc_mode[::1], fmt=f"ok", ecolor='black', elinewidth=0.5, ms=3)
            elif errbars == False:
                plt.errorbar(fs_mode[::1], vs_mode[::1], fmt=f"ok", elinewidth=0.5, ms=3)
        if display_method == "disp":
            plt.show()
        elif display_method == "save":
            plt.savefig(path, format='png', dpi='figure', bbox_inches='tight')
            plt.close()


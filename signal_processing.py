"""
Author : José CUNHA TEIXEIRA
License : SNCF Réseau, UMR 7619 METIS
Date : January 10, 2022
"""


import sys
import colormaps as cmaps
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector
from scipy.fft import rfft, rfftfreq
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

from display import diag_print

np.set_printoptions(threshold=sys.maxsize)




#----------------------------------------------------------------------------------------------------
def normalize(array, dt, clip_factor=6, clip_weight=10, norm_win=10, norm_method="1bit"):
    """
    Temporal normalization.
    """
    array2 = np.copy(array)
    for i, tr in enumerate(array2.T):
        if norm_method == 'clipping':
            lim = clip_factor * np.std(tr.data)
            tr[tr > lim] = lim
            tr[tr < -lim] = -lim
        elif norm_method == "clipping_iter":
            lim = clip_factor * np.std(np.abs(tr))
            # as long as still values left above the waterlevel, clip_weight
            while tr[np.abs(tr) > lim].size > 0:
                tr[tr > lim] /= clip_weight
                tr[tr < -lim] /= clip_weight
        elif norm_method == 'ramn':
            lwin = int(norm_win/dt)
            st = 0                                               # starting point
            N = lwin                                             # ending point
            while N < len(tr):
                win = tr[st:N]
                w = np.mean(np.abs(win)) / (2. * lwin + 1)
                # weight center of window
                tr[int(st + lwin / 2)] /= w
                # shift window
                st += 1
                N += 1
            # taper edges
            taper = get_window(len(tr))
            tr *= taper
        elif norm_method == "1bit":
            tr = np.sign(tr)
            tr = np.float32(tr)
    return array2




#----------------------------------------------------------------------------------------------------
def get_window(N, alpha=0.2):
    window = np.ones(N)
    x = np.linspace(-1., 1., N)
    ind1 = (abs(x) > 1 - alpha) * (x < 0)
    ind2 = (abs(x) > 1 - alpha) * (x > 0)
    window[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
    window[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
    return window



#----------------------------------------------------------------------------------------------------
def whiten(array, f_ech, freqmin, freqmax):

    array2 = np.empty(array.shape)
    nsamp = f_ech

    for i, tr in enumerate(array.T):
        
        n = len(tr)
        frange = float(freqmax) - float(freqmin)
        nsmo = int(np.fix(min(0.01, 0.5 * (frange)) * float(n) / nsamp))
        f = np.arange(n) * nsamp / (n - 1.)
        JJ = ((f > float(freqmin)) & (f<float(freqmax))).nonzero()[0]
            
        # signal FFT
        FFTs = np.fft.fft(tr)
        FFTsW = np.zeros(n) + 1j * np.zeros(n)

        # Apodization to the left with cos^2 (to smooth the discontinuities)
        smo1 = (np.cos(np.linspace(np.pi / 2, np.pi, nsmo+1))**2)
        FFTsW[JJ[0]:JJ[0]+nsmo+1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0]+nsmo+1]))

        # boxcar
        FFTsW[JJ[0]+nsmo+1:JJ[-1]-nsmo] = np.ones(len(JJ) - 2 * (nsmo+1))\
        * np.exp(1j * np.angle(FFTs[JJ[0]+nsmo+1:JJ[-1]-nsmo]))

        # Apodization to the right with cos^2 (to smooth the discontinuities)
        smo2 = (np.cos(np.linspace(0., np.pi/2., nsmo+1))**2.)
        espo = np.exp(1j * np.angle(FFTs[JJ[-1]-nsmo:JJ[-1]+1]))
        FFTsW[JJ[-1]-nsmo:JJ[-1]+1] = smo2 * espo

        whitedata = 2. * np.fft.ifft(FFTsW).real
        
        array2[:, i] = np.require(whitedata, dtype="float32")

    return array2





#----------------------------------------------------------------------------------------------------
def makeFV(XT, si, offset, vmin, vmax, dv, fmax):
    """   Constructs a FV dispersion diagram
    args :
        XT (numpy array) : data
        si (float) : sampling interval in seconds
        offsets (numpy array) : offsets in meter
        vmin, vmax (float) : velocities to scan in m/s
        dv (float) : velocity step in m/s
        fmax (float) : maximum frequency computed
    returns :
        fs : frequency axis
        vs : velocity axis
        FV: dispersion plot
    """
    Nx, Nt = XT.shape
    XF = rfft(XT, axis=(1), n=Nt)

    fs = rfftfreq(Nt, si)
    imax = np.where(fs > fmax)[0][0]
    fs = fs[0:imax]
    XF = XF[: , 0:imax]

    vs = np.arange(vmin, vmax+dv, dv)

    FV = np.zeros((len(vs), len(fs)))

    dphi = 2 * np.pi * offset[:,np.newaxis,np.newaxis] * fs[np.newaxis,np.newaxis,:] / vs[np.newaxis,:,np.newaxis]
    XF = XF[:,np.newaxis,:]
    FV = np.abs(np.sum(np.exp(1j*dphi)*XF/abs(XF), axis=0))

    return FV, vs, fs




#----------------------------------------------------------------------------------------------------
def pick(FV_arr, fs, vs, dx, Nx, smooth=True, norm_method=None):
    diag_print("Warning", "Dispersion curve picking", "Select the dispersion zone to pick")

    FV = np.copy(FV_arr)
    if norm_method == "Frequency":
        for i, col in enumerate(FV.T):
            FV[:, i] = FV[:, i] / np.max(col)
    elif norm_method == "Global":
        FV = FV / np.max(FV)

    def select_callback(x):
        global poly_coords
        poly_coords = x
        diag_print("Info", "Polygon coordiantes", str(x))

    def toggle_selector(event):
        if event.key == 't':
            for selector in selectors:
                name = type(selector).__name__
                if selector.active:
                    print(f'{name} deactivated.')
                    selector.set_active(False)
                else:
                    print(f'{name} activated.')
                    selector.set_active(True)
        elif event.key == 'enter':
            plt.close()
        else :
            print('No action key pressed.')

    selectors = []

    fig, ax = plt.subplots()
    fig.set_size_inches(16,9)
    fig.set_dpi(100)
    extent = [fs[0], fs[-1], vs[0], vs[-1]]
    plt.imshow(np.flipud(FV), cmap="turbo", extent=extent, aspect="auto")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Phase Velocity (m/s)")
    plt.title(f"Dispersion image - Selection of picking zone")
    plt.colorbar()
    plt.xticks(np.arange(min(fs), max(fs)+20.0, 20.0))
    v_ticks = np.arange(0, max(vs), 200.0)
    v_ticks[0] = 1
    plt.yticks(v_ticks)

    selectors.append(PolygonSelector(
        ax, select_callback, useblit=True))
    fig.canvas.mpl_connect('key_press_event', toggle_selector)

    plt.show()

    f_picked, v_picked = extract_curve(FV, fs, vs, poly_coords, smooth)
    dc = lorentzian_error(v_picked, f_picked, dx, Nx, 0.5)

    return f_picked, v_picked, dc




#----------------------------------------------------------------------------------------------------
def extract_curve(FV, fs, vs, poly_coords, smooth):
    """
    Extracts f-v dispersion curve from f-v dispersion diagram by aiming maximums

    args :
        FV (2D numpy array) : dispersion diagram
        fs (1D numpy array) : frequency axis
        vs (1D numpy array) : velocity axis
        start (tuple of floats) : starting coordinates (f,v) values
        end (tuple of floats) : ending coordinates (f,v) values

    returns :
        curve (1D numpy array[velocity]) : f-v dispersion curve
    """

    df = fs[1] - fs[0]
    dv = vs[1] - vs[0]
    idx = np.zeros((len(poly_coords), 2), dtype=int)
    for i, (f,v) in enumerate(poly_coords):
        idx[i][0] = int(v/dv)
        idx[i][1] = int(f/df)

    poly_path = Path(idx)
    x,y = np.mgrid[:FV.shape[0], :FV.shape[1]]
    coors = np.hstack((x.reshape(-1, 1), y.reshape(-1,1)))

    mask = poly_path.contains_points(coors)
    mask = mask.reshape(FV.shape)

    FV_masked = FV * mask

    f_picked = []
    v_picked =[]

    v_start_i = np.min(idx[:, 0])
    v_end_i = np.max(idx[:, 0])
    f_start_i = np.min(idx[:, 1])
    f_end_i = np.max(idx[:, 1])

    FV_tmp = FV_masked[v_start_i:v_end_i, f_start_i+1:f_end_i]

    for i, FV_f in enumerate(FV_tmp.T): # FV_f is a vector of velocities for a frequency f
        v_max_i = np.where(FV_f == FV_f.max())[0][0]
        v_max = vs[v_max_i+v_start_i]
        if v_max_i+v_start_i == v_end_i-1 and i != 0:
            v_picked.append(v_picked[-1])
        else:
            v_picked.append(v_max)
        f_picked.append(fs[i+f_start_i+1])

    f_picked = np.array(f_picked)
    v_picked = np.array(v_picked)

    if smooth == True:
        if (len(v_picked)/2) % 2 == 0:
            wl = len(v_picked)/2 + 1
        else:
            wl = len(v_picked)/2
        v_picked_curve = savgol_filter(v_picked, window_length=wl, polyorder=4, mode="nearest")
    else :
        v_picked_curve = v_picked

    plt.figure(figsize=(16, 9), dpi=100)
    extent = [fs[0], fs[-1], vs[0], vs[-1]]
    plt.imshow(np.flipud(FV_masked**2), cmap=cmaps.wh_bl_gr_ye_re, extent=extent, aspect="auto")
    plt.plot(f_picked, v_picked, '+w', f_picked, v_picked_curve, 'r')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Phase Velocity (m/s)")
    plt.title(f"Dispersion image with real and smoothed picked curve")
    plt.colorbar()
    plt.xticks(np.arange(min(fs), max(fs)+20.0, 20.0))
    v_ticks = np.arange(0, max(vs), 200.0)
    v_ticks[0] = 1
    plt.yticks(v_ticks)
    plt.show()

    return f_picked, v_picked_curve



#----------------------------------------------------------------------------------------------------
def resamp(f, v, err, type="wavelength"):
    w = v / f
    func_v = interp1d(w, v)
    func_err = interp1d(w, err)

    if type == "wavelength":
        w_resamp = np.linspace(min(w), max(w), len(f))
    elif type == "wavelength-log":
        w_resamp = np.geomspace(min(w), max(w), len(f))
    
    v_resamp = func_v(w_resamp)
    err_resamp = func_err(w_resamp)
    f_resamp = v_resamp/w_resamp
    
    return f_resamp[::-1], v_resamp[::-1], err_resamp[::-1]



#----------------------------------------------------------------------------------------------------
def lorentzian_error(v_picked, f_picked, dx, Nx, a):
    # Resolution
    Dc_left = 1 / (1/v_picked - 1/(2*f_picked*Nx*dx))
    Dc_right = 1 / (1/v_picked + 1/(2*f_picked*Nx*dx))
    Dc = np.abs(Dc_left - Dc_right)
    
    # Absolute uncertainty
    dc = (10**-a) * Dc

    for i, (err, v) in enumerate(zip(dc, v_picked)):
        if err > 0.6*v :
            dc[i] = 0.6*v

    return dc
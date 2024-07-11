"""
Author : José CUNHA TEIXEIRA
License : SNCF Réseau, UMR 7619 METIS
Date : November 30, 2023
"""

import sys
import numpy as np
from os import mkdir, path, system
from time import time
from obspy import read
from scipy.fft import fft, fftfreq, rfft, rfftfreq
from scipy.signal import butter, correlate, filtfilt
from scipy.signal.windows import tukey
from display import display_dispersion_img, display_spectrum_img_fromArray, display_seismic_wiggle_fromStream, diag_print
from obspy2numpy import array_to_stream, stream_to_array
from signal_processing import normalize, whiten, makeFV
# from stackmaster.core import pws, tfpws, tfpws_dost, robust, adaptive_filter, selective, clusterstack
from folders import results_dir, data_dir
from slant_stack import slant_stack


### FUNCTIONS -------------------------------------------------------------------------------------
def cut(TX_raw, ts, t_cut_start, t_cut2):
    i_cut_start = np.where(ts >= t_cut_start)[0][0]
    i_cut2 = np.where(ts >= t_cut2)[0][0]
    ts_cut = np.arange(0, t_cut2-t_cut_start, ts[1]-ts[0])
    TX = TX_raw[i_cut_start:i_cut2,::]
    for i, data in enumerate(np.transpose(TX)):
        TX[:,i] = np.transpose(data * tukey(len(ts_cut)))
    return TX, ts_cut
### -----------------------------------------------------------------------------------------------




### DATA DIRECTORIES ------------------------------------------------------------------------------
profile = "P1" #Profile number
data_dir = f"{data_dir}/" + profile + "/Passive/"
### -----------------------------------------------------------------------------------------------




### RECORDS PARAMS --------------------------------------------------------------------------------
if profile == 'P1' :
    data_files = ["101", "102", "103", "104", "105", "106"]
    durations = [130, 90, 90, 90, 90, 90]
    FK_ratio_min = 0.3

if profile == 'P2':
    data_files = ["201", "202", "203", "204", "205", "206", "207", "208", "209", "210"]
    durations = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
    FK_ratio_min = 0.6

if profile == 'P3':
    data_files = ["301", "302", "303", "304", "305", "306", "307", "308", "309", "310"]
    durations = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
    FK_ratio_min = 0.6

delta_cut = 5
window_step = 1

cuts = []
files = []
for file, duration in zip(data_files, durations) :
    cuts_list = list(np.arange(0,duration-delta_cut+window_step, window_step))
    N = len(cuts_list)
    cuts += cuts_list
    cuts[-1] = cuts[-1]-0.002
    files += [file]*N

# Use first record for initialisation
stream = read(data_dir + files[0] + '.dat')
### -----------------------------------------------------------------------------------------------




### SENSORS PARAMS --------------------------------------------------------------------------------
N_sensors = 96
sensor_start = 0
sensor_end = sensor_start + N_sensors - 1
if sensor_end - sensor_start != N_sensors - 1 :
    diag_print("Error", "passive_processing", "sensor_start - sensor_end != N_sensors - 1")
    raise SystemExit

stream = stream[sensor_start:sensor_end+1]
Nx = len(stream)
if Nx != N_sensors :
    diag_print("Error", "passive_processing", "Nx != N_sensors")
    raise SystemExit

x_sensors = np.arange(0.0, 24, 0.25)

x_sensors = x_sensors[sensor_start:sensor_end+1]
if len(x_sensors) != N_sensors :
    diag_print("Error", "passive_processing", "len(x_sensors) != N_sensors")
    raise SystemExit
### -----------------------------------------------------------------------------------------------




### INTERFEROMETRY PARAMS -------------------------------------------------------------------------
virtual_sources = [1, N_sensors]
### -----------------------------------------------------------------------------------------------




### MASW PARAMS -----------------------------------------------------------------------------------
W_MASW = N_sensors
start_MASW = sensor_start
### -----------------------------------------------------------------------------------------------




### DISPLAY PARAMS --------------------------------------------------------------------------------
display_method = "save"
### -----------------------------------------------------------------------------------------------




### INITIALISATION --------------------------------------------------------------------------------
fmin = 0.0
fmax = 200.0
fs = None
vmin = 1.0
vmax = 1000.0
vs = None
dt = 0.002
shot_length = 0.25
interf_length = 2.0
interf_db = np.zeros((len(virtual_sources), Nx, len(files), int(shot_length/dt)))
interf_db_stack = np.zeros((Nx, int(interf_length/dt)))
### -----------------------------------------------------------------------------------------------




### FK SELECTION ----------------------------------------------------------------------------------
N_segs_max = None
FK_ratios = []
### -----------------------------------------------------------------------------------------------




### PWS -------------------------------------------------------------------------------------------
pws_nu = 2
### -----------------------------------------------------------------------------------------------




### RESULTS DIRECTORIES ---------------------------------------------------------------------------
if not path.exists(results_dir):
    mkdir(results_dir)
results_dir = f"{results_dir}/" + profile
if not path.exists(results_dir):
    mkdir(results_dir)
results_dir = f"{results_dir}/" + "Passive"
if not path.exists(results_dir):
    mkdir(results_dir)
if N_segs_max == None :
    results_dir = f"{results_dir}/" + f"all-Sl{delta_cut}-Ss{window_step}-FK{int(FK_ratio_min*100)}-PWS{pws_nu}/"
elif N_segs_max != None :
    results_dir = f"{results_dir}/" + f"all-Sl{delta_cut}-Ss{window_step}-FK{int(FK_ratio_min*100)}-Nsegsmax{N_segs_max}-PWS{pws_nu}/"
if not path.exists(results_dir):
        mkdir(results_dir)

results_dir = results_dir + f"W{W_MASW}/"
if not path.exists(results_dir):
        mkdir(results_dir)

if not path.exists(results_dir):
    mkdir(results_dir)

file_path = results_dir + f"{profile}_FKratios.txt"
if path.exists(file_path):
    system(f"rm {file_path}")
file_FK_ratios = open(file_path, "a")
### -----------------------------------------------------------------------------------------------




tic = time()
for segment_i, file in enumerate(files):
    cut_start = cuts[segment_i]

    diag_print("Info", f"Main", f"Processing file {file} - cut {cut_start} - delta_cut {delta_cut}")

    stream = read(data_dir + file + ".dat")
    stream = stream[sensor_start:sensor_end+1]
    stream.detrend('demean')
    stream.detrend("linear")



    ### TIME PARAMS -------------------------------------------------------------------------------
    if stream[0].stats.delta != dt:
        ratio = int(dt/ stream[0].stats.delta)
        stream.decimate(ratio)
    Nt = stream[0].stats.npts
    dt = stream[0].stats.delta
    ts_raw = np.arange(0, Nt*dt, dt)
    f_ech = stream[0].stats.sampling_rate
    ### -------------------------------------------------------------------------------------------



    ### RAW DATA ----------------------------------------------------------------------------------
    # Data array ---
    TX_raw = stream_to_array(stream, len(stream), Nt)

    # Remove trace ---
    # TX_raw[:,16] = np.zeros((Nt))

    # Filter ---
    lowcut = 49
    highcut = 51
    order = 4
    nyq = 0.5 * f_ech
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandstop')
    for i, trace in enumerate(TX_raw.T):
        TX_raw[:,i] = filtfilt(b, a, trace)

    lowcut = 99
    highcut = 101
    order = 4
    nyq = 0.5 * f_ech
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandstop')
    for i, trace in enumerate(TX_raw.T):
        TX_raw[:,i] = filtfilt(b, a, trace)

    lowcut = 149
    highcut = 151
    order = 4
    nyq = 0.5 * f_ech
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandstop')
    for i, trace in enumerate(TX_raw.T):
        TX_raw[:,i] = filtfilt(b, a, trace)
    ### -------------------------------------------------------------------------------------------



    ### CUT ---------------------------------------------------------------------------------------
    TX, _ = cut(TX_raw, ts_raw, cut_start, cut_start+delta_cut)
    ### -------------------------------------------------------------------------------------------



    ### FK SELECTION (Cheng et al., 2018) ---------------------------------------------------------
    fs = rfftfreq(TX.shape[0], dt)
    fimax = np.where(fs > fmax)[0][0]
    fimin = np.where(fs > fmin)[0][0]
    fs = fs[fimin:fimax]

    ks = fftfreq(TX.shape[1], x_sensors[1]-x_sensors[0])

    FX = rfft(TX, axis=0)
    FX = FX[fimin:fimax,:]
    n = FX.shape[1]

    FK = fft(FX, axis=1)

    if n%2 == 0:
        K_minus = abs(FK[:, 1:int(n/2+1)])
        K_plus = abs(FK[:, int(n/2):])
    else:
        K_minus = abs(FK[:, 1:int((n-1)/2+1)])
        K_plus = abs(FK[:, int(((n-1)/2)+1):])

    FK_ratio = 0
    source_direction = None

    if np.sum(K_plus) > np.sum(K_minus):
        FK_ratio = np.sum(K_plus)/np.sum(K_minus) - 1
        if FK_ratio > FK_ratio_min :
            source_direction = "L"
            FK_ratios.append([segment_i, abs(FK_ratio)])

    elif np.sum(K_minus) > np.sum(K_plus):
        FK_ratio = - np.sum(K_minus)/np.sum(K_plus) + 1
        if FK_ratio < -FK_ratio_min :
            source_direction = "R"
            FK_ratios.append([segment_i, abs(FK_ratio)])

    print(source_direction, FK_ratio)

    file_FK_ratios.write(f"{FK_ratio}\n")
    ### -------------------------------------------------------------------------------------------


    if source_direction == "L" or source_direction == "R" :
        ### TEMPORAL NORMALIZATION (optional) -----------------------------------------------------
        # TX = normalize(TX, dt, clip_factor=1.0, norm_method="clipping")
        # TX = normalize(TX, dt, norm_win=0.01, norm_method="ramn")
        ### ---------------------------------------------------------------------------------------



        ### WHITENING -----------------------------------------------------------------------------
        TX = whiten(TX, f_ech, 0, fmax)
        ### ---------------------------------------------------------------------------------------



        ### INTERFEROMETRY ------------------------------------------------------------------------
        for i_s, virtual_source in enumerate(virtual_sources) :
            t_s = TX.T[virtual_source-1]

            for i_r, t_r in enumerate(TX.T):
                correl = correlate(t_r, t_s)
                tmp = np.copy(correl)
                tmp_0 = tmp[int(len(correl)//2)]
                tmp_del_tmp_0 = np.delete(tmp, int(len(correl)//2))
                acausal, causal = np.hsplit(tmp_del_tmp_0, 2)
                acausal = np.flip(acausal)

                if source_direction == "L" :
                    if virtual_source == 1:
                        correl_sym = causal
                    elif virtual_source == N_sensors:
                        correl_sym = acausal
                
                elif source_direction == "R" :
                    if virtual_source == 1:
                        correl_sym = acausal
                    elif virtual_source == N_sensors:
                        correl_sym = causal

                correl_sym = np.insert(correl_sym, 0, tmp_0)
                correl_sym = correl_sym[0:int(shot_length/dt)]
                
                for i, (correl_val, tukey_val)  in enumerate(zip(correl_sym, tukey(len(correl_sym)))):
                    if i >= len(correl_sym)//2:
                        correl_sym[i] = correl_val * tukey_val
                interf_db[i_s, i_r, segment_i, :] = correl_sym




### STACK INTERFEROMETRY --------------------------------------------------------------------------
diag_print("Info", f"Main", f"Interferometry stack")

FK_ratios = np.array(FK_ratios)
FK_ratios = FK_ratios[FK_ratios[:, 1].argsort()]

print(f"Interf database size : {interf_db.shape}")
# Remove empty virtual gathers (source_direction == None)
index = []
for i in range(len(files)) :
    if not np.any(interf_db[:, :, i, :]) :
        index.append(i)
if N_segs_max != None :
    if len(FK_ratios) > N_segs_max :
        for i in range(len(FK_ratios) - N_segs_max):
            index.append(int(FK_ratios[i,0]))
interf_db = np.delete(interf_db, index, 2)
print(f"Interf database size after removal of not selected files: {interf_db.shape}\n")


for i_r in range(Nx) :
    if len(virtual_sources) == 2 :
        arr1 = np.copy(interf_db[0, i_r, :, :])
        arr2 = np.copy(interf_db[1, Nx-1-i_r, :, :])
        arr = np.vstack((arr1, arr2))
    
    elif len(virtual_sources) == 1:
        if virtual_sources[0] == 1:
            arr = np.copy(interf_db[0, i_r, :, :])
        elif virtual_sources[0] == N_sensors:
            arr = np.copy(interf_db[0, Nx-1-i_r, :, :])

    index = []
    for i, tr in enumerate(arr):
        if not np.any(tr):
            index.append(i)
    arr = np.delete(arr, index, 0)

    # Linear stack ---
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = np.sum(arr, axis=0) / len(cuts)

    # PWS ---
    st = array_to_stream(arr.T, dt, range(arr.shape[0]))
    st = st.stack(stack_type=("pw", pws_nu))
    interf_db_stack[i_r, 0:int(shot_length/dt)] = st[0].data

    # Other stacks from stackmaster
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = tfpws(np.hstack((arr, np.zeros((arr.shape[0], int(1.75/dt))))),p=2)[0:int(shot_length/dt)]
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = tfpws(arr,p=5)
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = tfpws_dost(arr,p=2)
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = robust(arr,epsilon=1E-5,maxstep=10,win=None,stat=False,ref=None)
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = adaptive_filter(arr,g=1)
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = selective(arr,cc_min=0,epsilon=1E-5,maxstep=10,win=None,stat=False,ref=None)
    # interf_db_stack[i_r, 0:int(shot_length/dt)] = clusterstack(arr,h=0.90,win=None,axis=0,normalize=True,plot=False)

# Stream format ---
offsets = np.abs(x_sensors - x_sensors[0])
st_interf = array_to_stream(interf_db_stack.T, dt, offsets)
name_path = results_dir +  f"{profile}_VSG.segy"
st_interf.write(name_path, format="SEGY", data_encoding=1, byteorder=sys.byteorder)

name_path = results_dir + f"{profile}_VSG.png"
display_seismic_wiggle_fromStream(st_interf, x_sensors, display_method=display_method, path=name_path, norm_method="trace")

# Spectrum ---
name_path1 = results_dir + f"{profile}_spectrogram.png"
name_path2 = results_dir + f"{profile}_spectrogramFirstLastTrace.png"
display_spectrum_img_fromArray(interf_db_stack.T, dt, offsets, display_method=display_method, path1=name_path1, path2=name_path2, norm_method="trace")
### -----------------------------------------------------------------------------------------------




### SLANT STACK -----------------------------------------------------------------------------------
diag_print("Info", f"Main", f"Slant Stack")

arr = np.copy(interf_db_stack)

offsets = np.abs(x_sensors - x_sensors[0])
(FV, vs, fs) = slant_stack(arr, dt, offsets, vmin, vmax+0.1, 0.1, fmax)

name_path = results_dir + f"{profile}_dispIm.png"
display_dispersion_img(FV, fs, vs, display_method='save', path=name_path, normalization='Frequency')

name_path = results_dir + f"{profile}_dispIm.npy"
np.save(name_path, FV)

name_path = results_dir + f"{profile}_fs.npy"
np.save(name_path, fs)

name_path = results_dir + f"{profile}_vs.npy"
np.save(name_path, vs)
### -----------------------------------------------------------------------------------------------



### END -------------------------------------------------------------------------------------------
file_FK_ratios.close()

file_log = open(results_dir + f"{profile}_log.txt", "a")
file_log.write(f"Profil : {profile}\n")
file_log.write(f"Data directory : {data_dir}\n\n")
file_log.write(f"Data files : {data_files}\n")
file_log.write(f"Durations : {durations}\n")
file_log.write(f"Segment length : {delta_cut}\n")
file_log.write(f"Segment step : {window_step}\n\n")
file_log.write(f"Number of sensors : {N_sensors} | from sensor {sensor_start} to sensor {sensor_end}\n")
file_log.write(f"Sensor positions (m) : {x_sensors}\n\n")
file_log.write(f"Virtual sources : {virtual_sources}\n\n")
file_log.write(f"MASW window : {W_MASW} | MASW start : {start_MASW}\n")
file_log.write(f"fmin : {fmin} | fmax : {fmax} | fs : {fs}\n")
file_log.write(f"vmin : {vmin} | vmax : {vmax} | vs : {vs}\n\n")
file_log.write(f"Sample spacing (s) : {dt} | shot length (s) : {shot_length} | interf length : {interf_length}\n\n")
file_log.write(f"Minimum FK ratio : {FK_ratio_min} | Maximum segment number : {N_segs_max} | Number of used segmnts : {interf_db.shape[2]}\n\n")
file_log.write(f"PWS power : {pws_nu}\n\n")
file_log.close()

toc = time()
print(f"Elapsed Time : {toc-tic} s")
### -----------------------------------------------------------------------------------------------

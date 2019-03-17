import numpy as np
from lpf import lowpass_cosine as filt
import matplotlib.pyplot as plt
plt.ion()

def fwhm(v_array):

    tau = 1.0
    fc = 0.01
    info = np.loadtxt("./HI_INFO", dtype = "str")
    Vmin = np.float(info[np.where(info == 'VMIN')[0][0]][1])
    Vmax = np.float(info[np.where(info == 'VMAX')[0][0]][1])
    steps = v_array.shape[1]
    Vres = (Vmax - Vmin)/steps
    fwhm = np.zeros(v_array.shape[0])
    vlpf = np.zeros((v_array.shape[0], v_array.shape[1] - 1))
    #near = np.zeros_like(fwhm).astype('int')
    #second = np.zeros_like(fwhm).astype('int')

    for i in range(len(v_array)):
        if not np.nanmean(v_array[i]):
            fwhm[i] = np.nan
        else:
            try:
                v_lpf = filt(v_array[i], tau, fc, fc/3.)
                vlpf[i] = v_lpf
                avg_v = np.mean(v_lpf)
                idxs = np.where(v_lpf <= avg_v)[0]
                d_idxs = np.diff(idxs)
                left_idx = idxs[np.argmax(d_idxs)]
                right_idx = idxs[left_idx + 1]
                width = right_idx - left_idx
                dV = width * Vres
                disp = np.sqrt(8.0*np.log(2))*dV
                dV /= 1.0e3
                disp /= 1.0e3
                print dV, disp
                #print "FWHM:", fwhm[i]
                #print
            except ValueError:
                fwhm[i] = np.nan
    return vlpf, fwhm

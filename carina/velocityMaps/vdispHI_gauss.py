import numpy as np
from lpf import lowpass_cosine as filt
import matplotlib.pyplot as plt
from scipy.stats import norm
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
    #near = np.zeros_like(fwhm).astype('int')
    #second = np.zeros_like(fwhm).astype('int')

    for i in range(len(v_array)):
        if not np.nanmean(v_array[i]):
            fwhm[i] = np.nan
        else:
            try:
                #v_lpf = filt(v_array[i], tau, fc, fc/3.)
                 
                print dV, disp
                #print "FWHM:", fwhm[i]
                #print
            except ValueError:
                fwhm[i] = np.nan
    return vlpf, fwhm

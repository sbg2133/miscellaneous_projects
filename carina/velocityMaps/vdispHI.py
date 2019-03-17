import numpy as np
from lpf import lowpass_cosine as filt
import matplotlib.pyplot as plt
plt.ion()

def fwhm(v_array, vmin, vmax):
    steps = v_array.shape[1]
    vres = (vmax - vmin)/steps

    v = np.arange(vmin, vmax, vres)/1.0e3

    tau = 1.0
    fc = 0.03

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
                maximum = np.nanmax(v_lpf)
                max_idx = np.nanargmax(v_lpf)
                minimum = np.nanmin(v_lpf)
                height = maximum - minimum
                half_max = height/2.
                nearest = np.nanargmin((np.abs(v_lpf - half_max)))
                #near[i] = nearest
                print "half max:", half_max
                print "nearest:", nearest
                print "max idx:", max_idx
                if nearest > max_idx: # if half max point is to right of maximum
                    second_point = np.nanargmin((np.abs(v_lpf[:max_idx] - half_max)))
                if nearest < max_idx: # if it's on the left
                    second_point = np.nanargmin((np.abs(v_lpf[max_idx:] - half_max)))
                #second[i] = second_point
                #print "second:", second_point
                half_width_idxs = np.abs(nearest - second_point)
                if half_width_idxs > 1000:
                    fwhm[i] = np.nan
                else:
                    fwhm[i] = 2.*half_width_idxs*vres
                print "FWHM:", fwhm[i]
                print
            except ValueError:
                fwhm[i] = np.nan
    return vlpf, fwhm

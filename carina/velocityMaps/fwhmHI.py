import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from lpf import lowpass_cosine as filt

def vfwhm(v_array):
    info = np.loadtxt("./HI_INFO", dtype = "str")
    Vmin = np.float(info[np.where(info == 'VMIN')[0][0]][1])
    Vmax = np.float(info[np.where(info == 'VMAX')[0][0]][1])
    steps = v_array.shape[1]
    Vres = (Vmax - Vmin)/steps
    V = np.arange(Vmin, Vmax, Vres)/1.0e3
    tau = 1.0
    fc = 0.01
    V = V[:-1]
    zero_idx = (np.abs(V-0.0)).argmin()
    fwhm = np.zeros(len(v_array))
    for i in range(len(v_array)):
        data = filt(v_array[i], tau, fc, fc/3.)
        data = data[zero_idx - 100: zero_idx + 150]
        mean, sigma = norm.fit(data)
        sigma *= (Vres/1.0e3)
        fwhm[i] = np.sqrt(8.0*np.log(2.0))*sigma
    return fwhm

#np.save('v_fwhm_HI_5arcmin_kms.npy', fwhm)
#np.save('v_fwhm_HI_nskip2.npy', fwhm)
"""
data = filt(v_array[0], tau, fc, fc/3.)
v = V[zero_idx - 100: zero_idx + 150]
data = data[zero_idx - 100: zero_idx + 150]
plt.plot(v, data/np.sum(data))
mean, sigma = norm.fit(data)
y = norm.pdf(v, mean, sigma)
plt.plot(v, y)
plt.show()
"""

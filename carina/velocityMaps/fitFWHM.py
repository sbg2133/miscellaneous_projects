import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel
from lmfit.models import GaussianModel
from lpf import lowpass_cosine as filt

def vfwhm(v_array, lor = True):
    info = np.loadtxt("./HI_INFO", dtype = "str")
    Vmin = np.float(info[np.where(info == 'VMIN')[0][0]][1])
    Vmax = np.float(info[np.where(info == 'VMAX')[0][0]][1])
    steps = v_array.shape[1]
    Vres = (Vmax - Vmin)/steps
    V = np.arange(Vmin, Vmax, Vres)/1.0e3
    tau = 1.0
    fc = 0.01
    #V = V[:-1]
    zero_idx = (np.abs(V-0.0)).argmin()
    x = V[zero_idx - 100: zero_idx + 150]
    fwhm = np.zeros(len(v_array))
    count = 0
    for i in range(len(v_array)):
        print len(v_array) - count
        #data = filt(v_array[i], tau, fc, fc/3.)
        y = v_array[i][zero_idx - 100: zero_idx + 150]
        # if lor = False, fit Guassian
        if not lor:
           mod = GaussianModel()
        else:
           mod = LorentzianModel()
        pars = mod.guess(y, x=x)
        out = mod.fit(y, pars, x=x)
        sigma = out.best_values['sigma']
        sigma *= (Vres/1.0e3)
        if not lor:
            fwhm[i] = np.sqrt(8.0*np.log(2.0))*sigma
        else:
            fwhm[i] = 2.0*sigma
        #print fwhm[i]
        count += 1
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

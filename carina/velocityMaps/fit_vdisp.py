import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from lpf import lowpass_cosine as filt

info = np.loadtxt("./HI_INFO", dtype = "str")
Vmin = np.float(info[np.where(info == 'VMIN')[0][0]][1])
Vmax = np.float(info[np.where(info == 'VMAX')[0][0]][1])
v_array = np.load('./large_files/v_array_HI_5.npy')
steps = v_array.shape[1]
Vres = (Vmax - Vmin)/steps
V = np.arange(Vmin, Vmax, Vres)/1.0e3
tau = 1.0
fc = 0.02
V = V[:-1]
for i in range(len(v_array)):
    data = filt(v_array[i], tau, fc, fc/3.)
    mean, sigma = norm.fit(data)
    fwhm = np.sqrt(8.0*np.log(2.0))*sigma
    print sigma, fwhm

#data = filt(v_array[0], tau, fc, fc/3.)
plt.plot(data/np.sum(data))
mean, sigma = norm.fit(data)
y = norm.pdf(mean, sigma)
plt.plot(y)
plt.show()

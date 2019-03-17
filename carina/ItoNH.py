import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy
from astropy.wcs import WCS
import sys, os
from getIQU import IQU
from astropy import coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.interpolate import griddata
plt.ion()

root_dir = '/home/wizwit/miscellaneous_projects/carina/carinaData'
blast250_file = os.path.join(root_dir, 'smooth/3.0_arcmin/carinaneb_250_smoothed_3.0_rl.fits')

beta = 1.27

def getPsi(path_to_file):
    I, Q, U, __, wcs = IQU(path_to_file)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, psi

I, __, __, wcs_250, __, = getPsi(blast250_file)

#tau_d = (nu/nu0)**beta
# See Walker pg. 71
# nu0 = frequency at which dust emission becomes optically thin
#nu0 = 0.103 * Td # 0.103 (THz/K) * Td
#Inu_dust = Bnu(Td)*(1.0 - np.exp(1.0 - e**(-1.0*tau_d))

# See Walker pg. 69
# Av = 1.086*tau_d
# N_H = 1.79e21 * Av # (atoms/cm**2 mag)

# 1) Solve tau_d for temperature
# 2) Plug into Inu_dust equation





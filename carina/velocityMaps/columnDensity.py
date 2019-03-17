import numpy as np
import sys, os
from astropy.io import fits

h = 6.626e-27
c = 3.0e10
k = 1.38e-16
D = 2.3e3 * 3.086e18 # distance to Carina, pc to cm
f = 100.0 # molecular gas to dust ratio
mH2 = 3.32e-24 # g

# Planck function
def B(nu, T): # erg / s * cm^2 * Hz * str
  n = (np.exp(h*nu/(k*T)) - 1.0)**-1.0
  return (2*h*nu**3 / c**2)*n

# dust opacity, cm^2/g
def kappa(um):
    return 13.16*(um/250.)**-2

# column depth, cm
def cloud_depth(pc):
    return 3.086e18 * pc

# column density, cm^-2
def NH2(band, I):
    if band == '250':
        nu = c/250.0e-4
        kapp = kappa(250)
    if band == '350':
        nu = c/350.0e-4
        kapp = kappa(350)
    if band == '500':
        nu = c/500.0e-4
        kapp = kappa(500)

    root_dir = "/home/wizwit/miscellaneous_projects/carina/carinaData/planckData"
    T_file = os.path.join(root_dir, "planck_353_T.fits")
    Av_file = os.path.join(root_dir, "planck_353_AvR.fits")
    Beta_file = os.path.join(root_dir, "planck_353_beta.fits")

    Td = fits.open(T_file)[1].data
    Av = fits.open(Av_file)[1].data
    beta = fits.open(Beta_file)[1].data
    Omega_beam = np.pi * ((5.0/60.) * (np.pi/180.0))**2 # steradians

    Bnu = B(nu, Td) / 1.0e-23 # Jy/str
    Fnu = (I*1.0e6) * Omega_beam # Jy
    tau = Fnu / (Omega_beam * Bnu)
    N = tau*f/(kapp*mH2)
    return N, tau, Av

# nH2, cm^-3
def nH2(band, I, L):
    N, tau, Av = NH2(band, I)
    n = N / cloud_depth(L) # thickness in pc
    n[n < 0] = np.nan
    N = N[30:-30,260:-260]
    N = N[N > 0]
    N = N[~np.isnan(N)]
    Av = N / 9.4e20
    t = tau[30:-30,260:-260]
    t = t[~np.isnan(t)]
    print "< NH2 > =", np.nanmean(N)
    print "< tau > =", np.nanmean(t)
    print "< Av > =", np.nanmean(Av)
    print "< nH2> =", np.nanmean(n)
    print "min nH2 =", np.nanmin(n)
    print "max nH2 =", np.nanmax(n)
    return n

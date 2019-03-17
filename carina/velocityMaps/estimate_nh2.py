import numpy as np
import sys, os
from getIQU import IQU
from subprocess import call
import sys, os
import glob
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.visualization import AsinhStretch
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import scipy.ndimage as ndimage
from skimage import filters
import aplpy
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, PowerStretch)
from scipy import constants as const
from makePretty import pretty
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ion()

band = sys.argv[1] # '250', '350', or '500'
L = float(sys.argv[2]) # A number of pc (e.g., 0.5, 5)
print L

h = const.h
h = 6.26e-27
c = const.c
c = 3.0e10
k = const.k
k = 1.38e-16

D = 2.3e3 * 3.086e18 # distance to Carina, pc to cm

def cloud_depth(pc):
    return 3.086e18 * pc

#L = D*np.tan(np.deg2rad(0.5)) # estimated depth of cloud, ~0.5 deg
a = 1.0e-4 # grain radius, cm
rho = 3.0 # grain density, 3 g / cm^3

def kappa(um):
    return 13.16*(um/250.)**-2

rootdir = "/home/wizwit/miscellaneous_projects/carina/carinaData/smooth/5.0_arcmin_no_kernnorm"
if band == '250':
    blast_file = os.path.join(rootdir, "carinaneb_250_smoothed_5.0_rl_cal.fits")
    kapp = kappa(250)
if band == '350':
    blast_file = os.path.join(rootdir, "carinaneb_350_smoothed_5.0_rl_cal.fits")
    kapp = kappa(350)
if band == '500':
    kapp = kappa(500)
    blast_file = os.path.join(rootdir, "carinaneb_500_smoothed_5.0_rl_cal.fits")

I, __, __, __, wcs = IQU(blast_file)
cen_coord = [160.84429, -59.582361]

#I = I[30:-30,260:-260]

pix_deg = 0.00277777777778 # increment in degrees at center of image
xsize_deg = I.shape[0]*pix_deg
ysize_deg = I.shape[1]*pix_deg
# str
Omega_source = xsize_deg * ysize_deg * (np.pi/180.0)**2

root_dir = "/home/wizwit/miscellaneous_projects/carina/carinaData/planckData"
T_file = os.path.join(root_dir, "planck_353_T.fits")
#Av_file = os.path.join(root_dir, "planck_353_Av.fits")
Av_file = os.path.join(root_dir, "planck_353_AvR.fits")
Beta_file = os.path.join(root_dir, "planck_353_beta.fits")

#Td = fits.open(T_file)[1].data[30:-30,260:-260]
#Av = fits.open(Av_file)[1].data[30:-30,260:-260]
#beta = fits.open(Beta_file)[1].data[30:-30,260:-260]
Td = fits.open(T_file)[1].data
Av = fits.open(Av_file)[1].data
beta = fits.open(Beta_file)[1].data

#Omega_beam = np.pi* (3.0 * 60.0 / 206265.0)**2 # arcmin * 60 arcsec / arcmin * 1 rad / 206265 arcsec
Omega_beam = np.pi * ((5.0/60.) * (np.pi/180.0))**2
NH2 = 9.4e20/Av # 1/cm^2 See Franco 2015, Bohlin 1978

if band == '250':
    nu = c/250.0e-4
if band == '350':
    nu = c/350.0e-4
if band == '500':
    nu = c/500.0e-4

nu0 = 0.103 * Td * 1.0e12 # nu0 = 0.103 THz/K See Walker pg 71
Qv = 7.5e-4 * ((nu/1.0e12) / 2.4 )**beta # Walker, pg 74
tau_d = (nu/nu0)**beta
f = 100.0 # molecular gas to dust ratio

# Planck function
def B(nu, T): # erg / s * cm^2 * Hz * str
  n = (np.exp(h*nu/(k*T)) - 1.0)**-1.0
  return (2*h*nu**3 / c**2)*n

Bnu = B(nu, Td) / 1.0e-23 # Jy/str

Fnu = (I*1.0e6) * Omega_beam # Jy
tau = Fnu / (Omega_beam * Bnu)
tau_Av = Av/1.086

Msun = 1.989e32 # g
# Dust mass
Md = (4./3.)*((D**2)*a*rho*tau_Av*Omega_beam)/(Qv) # g
#Md = (4./3.)*((D**2)*a*rho*Fnu/B)/(Qv) # g
MH2 = f*Md

# NH2
mH2 = 3.32e-24 # g
#NH2_est = np.pi*MH2 / (4.0*mH2*(D**2.0)*Omega_beam**2.0)
NH2_est = tau*f/(kapp*mH2)
Av_est = NH2_est / 9.4e20

# nH2, kg / m^3
#nH2_est = 6.0 * f * Md / (np.pi * D**3.0 * Omega_beam**3.0 * mH2)
#nH2_est2 = (4./3.)*6.0*f*Fnu*a*rho / (np.pi*Qv*B_250*D*Omega_beam**3.0)
#nH2_cgs = nH2_est * 1.0e3 * (1.0e2)**3.0

# nH2, 1/cm^3
print L
nH2 = NH2_est / cloud_depth(L) # thickness in pc
nH2[nH2 < 0] = np.nan

print "< nH2> =", np.nanmean(nH2)
print "min nH2 =", np.nanmin(nH2)
print "max nH2 =", np.nanmax(nH2)
print
#nH2 = NH2 / (D*Omega_beam)

np.save('nh2_' + str(band) + '_' + str(L) + 'pc.npy', nH2)

"""
hdu = fits.PrimaryHDU(data=NH2_est, header=wcs.to_header())
f = aplpy.FITSFigure(hdu, figsize = (12, 12))
f.set_theme('publication')
ax = plt.gca()
f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
ax.set_facecolor("k")
f.show_colorscale(cmap = 'viridis', vmin = 3.0e20, vmax = 6.0e21)
f.add_scalebar(15/60.) # deg
f.scalebar.set_color('white')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('white')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.show()
f.colorbar.set_axis_label_text(r'$N(H_{2})$ cm$^{-2}$')
f.colorbar.set_axis_label_font(size = 16)
f.frame.set_linewidth(1)  # points
f.frame.set_color('black')
"""

N = NH2_est[30:-30,260:-260]
N = N[N > 0]
N = N[~np.isnan(N)]

t = tau[30:-30,260:-260]
t = t[~np.isnan(t)]

A = 9.4e20 / N

print "< NH2 > =", np.nanmean(N)
print "< tau > =", np.nanmean(t)
print "< Av > =", np.nanmean(A)

"""
plt.hist(np.log10(N), bins = 1000, color = 'C0')
plt.xlabel(r"log$_{10}$N(H$_{2})$ cm$^{-2}$")
plt.ylabel('N')
plt.xlim(19.8, 22.5)
pretty()
plt.savefig(os.path.join(save_files_here, 'NH2_hist.eps'), format='eps', bbox_inches = 'tight')
"""

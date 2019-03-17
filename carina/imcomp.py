from getIQU import IQU
from subprocess import call
import sys
import numpy as np
import glob
import matplotlib.pyplot as plt
import aplpy
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
import scipy.ndimage
from skimage import filters
plt.ion()

if len(sys.argv) < 2:
    print "You must supply a band, 'e.g., carina_lic.py 250'"
    sys.exit()
else:
    band = sys.argv[1]

bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']

# define file paths
blastpol_dir = './carinaData/smooth/3.0_arcmin'
filename = glob.glob(blastpol_dir + '/carinaneb_' + band + '_smoothed_3.0_rl.fits')[0]

# load in I, Q, U for desired band
Ivals, Qvals, Uvals, wcs = IQU(band, filename)

I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

# Correct pvals as in Jamil's thesis, 5.7
#pvals[pvals > 0.5] = np.nan
pvals /= pol_eff[bands.index(band)]
phi = 0.5*np.arctan2(U,Q)
dx = np.cos(phi)
dy = np.sin(phi)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
xs, ys = np.meshgrid(X,Y)


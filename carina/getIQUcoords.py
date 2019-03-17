import sys
import numpy as np
import glob
import aplpy
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.wcs import WCS
import scipy.ndimage
from skimage import filter

# lists over which to iterate:
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']

# define file paths
blastpol_dir = './carinaData/blastData/'
filenames = glob.glob(blastpol_dir+'*p10_good_C_gls_map_cal.fits')

fwhm_orig = 36. # 250um fwhm, "
omega_orig = 1.13*fwhm_orig**2
smooth_fwhm = 4.5 * 60. # 4.5' in "
fwhm_norm = smooth_fwhm / fwhm_orig
sigma = fwhm_norm/2.3548
gauss_kernel = Gaussian2DKernel(sigma)
percent= 0.075 # empirically-tweaked

def IQU(band, smooth = True):
    for b, filename in enumerate(filenames):
        # Read in the HDUlist from the FITS file:
        if b == bands.index(band):
            print(filename)
            hdulist = fits.open(filename)
            # Print some info on what's in it:
            hdulist.info()
            for s, param in enumerate(stokes):
                if (stokes[s] == 'I'):
                    Ivals_raw = hdulist[s+1].data
                    wcs = WCS(hdulist[s+1].header)
                    minval = np.nanmin(Ivals_raw)
                    maxval = np.nanmax(Ivals_raw)
                    # Subtracting mean flux from dark region
                    far_region = Ivals_raw[120:350, 730:800]
                    Ivals = Ivals_raw - np.mean(far_region)
                    #Ivals[Ivals < np.std(far_region)] = np.nan
                    if smooth:
                        Ivals = convolve(Ivals, gauss_kernel)
                    Ivals[Ivals == 0.0] = np.nan
                if (stokes[s] == 'Q'):
                    Qvals_raw = hdulist[s+1].data
                    wcs = WCS(hdulist[s+1].header)
                    far_region = Qvals_raw[120:350, 730:800]
                    Qvals = Qvals_raw - np.mean(far_region)
                    #Qvals[Qvals < np.std(far_region)] = np.nan
                    if smooth:
                        Qvals = convolve(Qvals, gauss_kernel)
                    Qvals[Qvals == 0.0] = np.nan
                if (stokes[s] == 'U'):
                    Uvals_raw = hdulist[s+1].data
                    wcs = WCS(hdulist[s+1].header)
                    far_region = Uvals_raw[120:350, 730:800]
                    Uvals = Uvals_raw - np.mean(far_region)
                    #Uvals[Uvals < np.std(far_region)] = np.nan
                    if smooth:
                        Uvals = convolve(Uvals, gauss_kernel)
                    Uvals[Uvals == 0.0] = np.nan
    return Ivals, Qvals, Uvals, wcs

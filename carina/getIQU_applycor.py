import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

stokes = ['I', 'Q', 'U']

def IQU(filename, do_cov = False):
    print(filename)
    covar_file = filename[:-5] + '_cov.fits'
    hdulist = fits.open(filename)
    print hdulist.info()
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
            Ivals[Ivals == 0.0] = np.nan
        if (stokes[s] == 'Q'):
            Qvals_raw = hdulist[s+1].data
            wcs = WCS(hdulist[s+1].header)
            far_region = Qvals_raw[120:350, 730:800]
            Qvals = Qvals_raw - np.mean(far_region)
            #Qvals[Qvals < np.std(far_region)] = np.nan
            Qvals[Qvals == 0.0] = np.nan
        if (stokes[s] == 'U'):
            Uvals_raw = hdulist[s+1].data
            wcs = WCS(hdulist[s+1].header)
            far_region = Uvals_raw[120:350, 730:800]
            Uvals = Uvals_raw - np.mean(far_region)
            #Uvals[Uvals < np.std(far_region)] = np.nan
            Uvals[Uvals == 0.0] = np.nan
    if do_cov:
        covar_hdu = fits.open(covar_file)
        covar = covar_hdu[1].data
    else:
        covar = None
    return Ivals, Qvals, Uvals, covar, wcs

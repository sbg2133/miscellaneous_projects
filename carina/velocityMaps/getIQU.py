import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import blast.util

stokes = ['I', 'Q', 'U']

def IQU(filename, do_cov = False):
    print filename
    covar_file = filename[:-5] + '_cov_cal.fits'
    print covar_file
    band = filename[94:97]
    print band
    if band == '250':
        pol_eff = 0.81
    if band == '350':
        pol_eff = 0.79
    if band == '500':
        pol_eff = 0.82
    print covar_file
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
        new_map_hdu = fits.HDUList([hdulist[0]])
        I_hdu = hdulist[1].copy()
        I_hdu.data = Ivals
        new_map_hdu.append(I_hdu)
        Q_hdu = hdulist[1].copy()
        Q_hdu.data = Qvals
        new_map_hdu.append(Q_hdu)
        U_hdu = hdulist[1].copy()
        U_hdu.data = Uvals
        new_map_hdu.append(U_hdu)
        covar_hdu = fits.open(covar_file)
        new_hdu = blast.util.calc_pol_p_phi(new_map_hdu, cov_hdus=covar_hdu, pol_eff=pol_eff, var_norm=0.7368)
        P = new_hdu[1].data
        phi = new_hdu[2].data
        p = new_hdu[3].data/100.0
        sig_P = new_hdu[4].data
        sig_phi = new_hdu[5].data
        sig_p = new_hdu[6].data/100.0
        pol_data = np.array([P, phi, p, sig_P, sig_phi, sig_p])
    else:
        pol_data = None
    return Ivals, Qvals, Uvals, pol_data, wcs

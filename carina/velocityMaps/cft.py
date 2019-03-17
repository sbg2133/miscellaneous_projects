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
plt.ion()

rootdir = '/home/wizwit/miscellaneous_projects/carina/carinaData'
#cube_filename = os.path.join(rootdir, 'mopraData/G287_288.-2.0_0.5.12CO.fits')
sigma_file = os.path.join(rootdir,\
         'mopraData/G287_288.-2.0_0.5.12CO_sigma.fits')
blast250_file = os.path.join(rootdir, 'smooth/3.0_arcmin/carinaneb_250_smoothed_3.0_rl.fits')

#hdulist = fits.open(cube_filename)
#data_cube = hdulist[0].data
#wcs = WCS(hdulist[0].header)
#disp = np.nanstd(data_cube, axis = 0)

hdulist_co12 = fits.open(sigma_file)
sigma = hdulist_co12[0].data
wcs_co12 = WCS(hdulist_co12[0].header)

#hdulist_250 = fits.open(blast250_filename)
#blast250 = hdulist_250[0].data
#wcs_250 = WCS(hdulist_250[0].header)

def getPsi(path_to_file):
    I, Q, U, __, wcs = IQU(path_to_file)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, psi

I, __, __, wcs_250, psi = getPsi(blast250_file)

cdelta = 0.00277777777778 * 60 # deg * 60 arcmin/deg

nskip = np.int(np.round(1.0/cdelta)) # points per 3 arcmin
r = nskip*3

Ix = np.arange(I.shape[1])
Iy = np.arange(I.shape[0])
IX, IY = np.meshgrid(Ix, Iy)
ra, dec = wcs_250.all_pix2world(IX, IY, 0)
radec = SkyCoord(ra*u.degree, dec*u.degree, frame='fk5')

x_sig = np.arange(sigma.shape[1])
y_sig = np.arange(sigma.shape[0])
XSIG, YSIG = np.meshgrid(x_sig, y_sig)
ra_sig, dec_sig = wcs_co12.all_pix2world(XSIG, YSIG, 0)
radec_sig = SkyCoord(ra_sig*u.degree, dec_sig*u.degree, frame='galactic').fk5

blast_ra = np.asarray(radec.ra)
blast_dec = np.asarray(radec.dec)
x = np.arange(blast_ra.shape[0])[::nskip]
y = np.arange(blast_dec.shape[1])[::nskip]
X,Y = np.meshgrid(x,y)
ra_points, dec_points = wcs_250.all_pix2world(X, Y, 0)
radec_points = SkyCoord(ra_points*u.degree, dec_points*u.degree, frame='fk5')

#pol_disp = np.load('pol_disp.npy')
pol_disp = np.empty_like(blast_ra)
pol_disp[:] = np.nan

skip = (slice(50, -50, nskip), slice(50, -50, nskip))

for i in range(radec.shape[0]):
    for j in range(radec.shape[1]):
        mask_inner = (blast_ra - blast_ra[i,j])**2 + (blast_dec - blast_dec[i,j])**2 < (r)**2
        mask_outer = (blast_ra - blast_ra[i,j])**2 + (blast_dec - blast_dec[i,j])**2 <= (r)**2
        mask = mask_inner ^ mask_outer
        comp_points = psi[mask]
        center_point = psi[i,j]
        S2 = (center_point - comp_points)**2
        S2_debias = S2 - np.var(S2)
        if np.size(S2_debias):
            S = np.sqrt((1.0/len(S2_debias))*np.sum(S2_debias))
            pol_disp[x[i]][y[j]] = S
        else:
            pass

#hdu_250 = fits.PrimaryHDU(I, header=wcs_250.to_header())
#hdu_sigma = fits.PrimaryHDU(sigma, header=wcs_co12.to_header())
#fig = plt.figure()
#f_sigma = aplpy.FITSFigure(hdu_sigma, figure = fig)
#f_sigma.show_colorscale(cmap = 'inferno')
#f_250 = aplpy.FITSFigure(hdu_250, figure = fig)
#f_250.show_colorscale(cmap = 'inferno')
plt.figure()
plt.xlim(np.max(blast_ra), np.min(blast_ra))
plt.pcolormesh(blast_ra, blast_dec, I, cmap = "inferno")
plt.scatter(blast_ra[skip], blast_dec[skip], c = 'r', alpha = 0.5)
#plt.pcolormesh(blast_ra, blast_dec, pol_disp, alpha = 0.5)
#plt.pcolormesh(radec_sig.ra, radec_sig.dec, np.fliplr(sigma), cmap = "viridis", alpha = 0.5)


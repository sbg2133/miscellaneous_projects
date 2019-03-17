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
blast250_file = os.path.join(rootdir,\
         'smooth/3.0_arcmin/carinaneb_250_smoothed_3.0_rl.fits')

hdulist = fits.open(cube_filename)
data_cube = hdulist[0].data
wcs = WCS(hdulist[0].header)
disp = np.nanstd(data_cube, axis = 0)

hdulist_co12 = fits.open(sigma_file)
sigma = hdulist_co12[0].data
wcs_co12 = WCS(hdulist_co12[0].header)

sigma = disp
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

# points to use for S calculation
x = np.arange(psi.shape[0])[::nskip]
x = x[4:-4] # don't use edges
y = np.arange(psi.shape[1])[::nskip]
y = y[4:-4] # don't use edges
X,Y = np.meshgrid(x,y)

Ix = np.arange(I.shape[1])
Iy = np.arange(I.shape[0])
IX, IY = np.meshgrid(Ix, Iy)
ra, dec = wcs_250.all_pix2world(IX, IY, 0)
radec = SkyCoord(ra*u.degree, dec*u.degree, frame='fk5')
galcoords = radec.galactic

x_sig = np.arange(sigma.shape[1])
y_sig = np.arange(sigma.shape[0])
XSIG, YSIG = np.meshgrid(x_sig, y_sig)
ra_sig, dec_sig = wcs_co12.all_pix2world(XSIG, YSIG, 0)
radec_sig = SkyCoord(ra_sig*u.degree, dec_sig*u.degree, frame='galactic').fk5

# points to use for S calculation
x = np.arange(psi.shape[0])[::nskip]
x = x[4:-4] # don't use edges
y = np.arange(psi.shape[1])[::nskip]
y = y[4:-4] # don't use edges
X,Y = np.meshgrid(x,y)

#plt.imshow(I, origin = "lower", cmap = "viridis")
#plt.pcolormesh(Iy,Ix,I)
#plt.scatter(Y, X, c = 'r')

#pol_disp = np.load('pol_disp.npy')
blast_ra = np.asarray(radec.ra)
blast_dec = np.asarray(radec.dec)

pol_disp = np.empty(np.shape(radec))
pol_disp[:] = np.nan
r = 3.0/60 # 3 arcmin in deg

for i in range(len(x)):
    for j in range(len(y)):
        mask_inner = (blast_dec[np.newaxis,:]-y[j])**2 + (blast_ra[:,np.newaxis]-x[i])**2 < (r)**2
        mask_outer = (blast_dec[np.newaxis,:]-y[j])**2 + (blast_ra[:,np.newaxis]-x[i])**2 <= (r)**2
        mask = mask_inner ^ mask_outer
        comp_points = psi[mask]
        center_point = psi[x[i]][y[j]]
        S2 = (center_point - comp_points)**2
        S2_debias = S2 - np.var(S2)
        S = np.sqrt((1.0/len(S2_debias))*np.sum(S2_debias))
        pol_disp[x[i]][y[j]] = S

#hdu_250 = fits.PrimaryHDU(I, header=wcs_250.to_header())
#hdu_sigma = fits.PrimaryHDU(sigma, header=wcs_co12.to_header())
#fig = plt.figure()
#f_sigma = aplpy.FITSFigure(hdu_sigma, figure = fig)
#f_sigma.show_colorscale(cmap = 'inferno')
#f_250 = aplpy.FITSFigure(hdu_250, figure = fig)
#f_250.show_colorscale(cmap = 'inferno')
plt.figure()
plt.xlim(np.max(np.asarray(radec.ra)), np.min(np.asarray(radec.ra)))
plt.pcolormesh(radec.ra, radec.dec, I, cmap = "inferno")
plt.pcolormesh(radec.ra, radec.dec, pol_disp, alpha = 0.5)
plt.pcolormesh(radec_sig.ra, radec_sig.dec, np.fliplr(sigma), cmap = "viridis", alpha = 0.5)


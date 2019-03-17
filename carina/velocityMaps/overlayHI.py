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

info = np.loadtxt("./HI_INFO", dtype = "str")
root_dir = info[np.where(info == 'ROOT_DIR')[0][0]][1]
cube_file = os.path.join(root_dir,\
           info[np.where(info == 'CUBE_FILE')[0][0]][1])
blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST250_FILE')[0][0]][1])

vmin = float(info[np.where(info == 'VMIN')[0][0]][1])
vmax = float(info[np.where(info == 'VMAX')[0][0]][1])

cdelta_v = float(info[np.where(info == 'CDELTA_V')[0][0]][1])
cdelta_blast = float(info[np.where(info == 'CDELTA_BLAST')[0][0]][1])

pol_disp = np.load('./large_files/pol_disp.npy')
#vdisp = np.load('v_array_co13.npy')
#v_fwhm = np.load('fwhmCO13.npy')

hdulist_v = fits.open(cube_file)
data_cube = hdulist_v[0].data
dummy_v = hdulist_v[0].data[0]

wcs_v = WCS(naxis=2)

#wcs_v.naxis = np.array([1176, 1451])
wcs_v.wcs.crpix = [563.0, 665.0]
wcs_v.wcs.cdelt = np.array([-0.00277777798786, 0.00277777798786])
wcs_v.wcs.crval = [287.41231676, -1.14784038739]
wcs_v.wcs.ctype = ["GLON-SIN", "GLAT-SIN"]

def getPsi(path_to_file):
    I, Q, U, __, wcs = IQU(path_to_file)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, psi

I, __, __, wcs, psi = getPsi(blast_file)

cdelta_blast *= 60 # deg * 60 arcmin/deg
nskip_blast = np.int(np.round(1.0/cdelta_blast)) # points per arcmin
r_blast = nskip_blast

cdelta_v *= 60 # deg * 60 arcmin/deg
nskip_v = np.int(np.round(1.0/cdelta_v)) # points per arcmin
r_v = nskip_v

# points to use for S calculation
x = np.arange(psi.shape[1])[::nskip_blast]
x = x[30:-30] # don't use edges
y = np.arange(psi.shape[0])[::nskip_blast]
y = y[5:-5] # don't use edges
X,Y = np.meshgrid(x,y)
mask_ra, mask_dec = wcs.all_pix2world(X, Y, 0)
mask_gal_coords = SkyCoord(mask_ra*u.degree, mask_dec*u.degree, frame='fk5').galactic

Ix = np.arange(I.shape[1])
Iy = np.arange(I.shape[0])
IX, IY = np.meshgrid(Ix, Iy)
ra, dec = wcs.all_pix2world(IX, IY, 0)
radec = SkyCoord(ra*u.degree, dec*u.degree, frame='fk5')

x_v = np.arange(dummy_v.shape[1])
y_v = np.arange(dummy_v.shape[0])
XV, YV = np.meshgrid(x_v, y_v)
v_ra, v_dec = wcs_v.all_pix2world(XV, YV, 0)
radec_v = SkyCoord(v_ra*u.degree, v_dec*u.degree, frame='galactic').fk5

#pol_disp = np.empty((len(Ix), len(Iy)))
#pol_disp[:] = np.nan

"""
pol_coords = np.column_stack((X.ravel(), Y.ravel()))
for coord in pol_coords:
    center = psi[coord[1], coord[0]]
    mask_inner = (Ix[np.newaxis,:] - coord[0])**2 + (Iy[:,np.newaxis] - coord[1])**2 < (r_blast)**2
    mask_outer = (Ix[np.newaxis,:] - coord[0])**2 + (Iy[:,np.newaxis] - coord[1])**2 <= (r_blast)**2
    mask = mask_inner ^ mask_outer
    comp_points = psi[mask]
    center = psi[coord[1], coord[0]]
    S2 = (center - comp_points)**2
    S2_debias = S2 - np.var(S2)
    #S = np.sqrt((1.0/len(S2))*np.sum(S2))
    S = np.sqrt((1.0/len(S2_debias))*np.sum(S2_debias))
    #I[coord[1], coord[0]] = 1.0e3
    pol_disp[coord[0], coord[1]] = S
"""
v_x, v_y = wcs_v.all_world2pix(mask_gal_coords.l, mask_gal_coords.b, 0)
# integer pixel values
v_x = np.round(v_x).astype('int')
v_y = np.round(v_y).astype('int')

v_mask_ra, v_mask_dec = wcs_v.all_pix2world(v_x, v_y, 0)
v_mask_fk5 = SkyCoord(v_mask_ra*u.degree, v_mask_dec*u.degree, frame='galactic').fk5

# make velocity dispersion map
vdisp = np.empty_like(v_x)
vdisp[:] = np.nan

mask_coords = np.column_stack((v_x.ravel(), v_y.ravel()))
for coord in mask_coords[::20]:
#    if coord[0] > 120 or coord[1] < 0:
#        continue
    center = dummy_v[coord[1], coord[0]]
    mask = (x_v[np.newaxis,:] - coord[0])**2 + (y_v[:,np.newaxis] - coord[1])**2 <= (r_v)**2
    comp_points = dummy_v[mask]
    dummy_v[mask] = 1.0e8

plt.xlim(np.max(np.asarray(radec.ra)), np.min(np.asarray(radec.ra)))
plt.pcolormesh(radec.ra, radec.dec, I, cmap = "inferno")
#plt.pcolormesh(radec.ra, radec.dec, pol_disp, alpha = 0.5)
plt.pcolormesh(radec_v.ra, radec_v.dec, dummy_v, cmap = "viridis", alpha = 0.1)
#plt.scatter(co12_mask_fk5.ra, co12_mask_fk5.dec, c = 'r', alpha = 0.5)
#plt.scatter(mask_ra, mask_dec, c = 'g', alpha = 0.5)

#plt.pcolormesh(XCO12, YCO12, co12, cmap = "viridis", alpha = 0.5)
#plt.scatter(co12_x, co12_y, c = 'r', alpha = 0.5)

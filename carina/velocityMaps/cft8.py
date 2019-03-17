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
from vdisp import fwhm as fwhm
from makePretty import pretty
plt.ion()

rootdir = '/home/wizwit/miscellaneous_projects/carina/carinaData'
cube_filename = os.path.join(rootdir, 'smooth/co12_smooth_3.0_arcmin.fits')
sig_file = os.path.join(rootdir,\
         'mopraData/G287_288.-2.0_0.5.12CO_sigma.fits')
blast250_file = os.path.join(rootdir,\
         'smooth/3.0_arcmin/carinaneb_250_smoothed_3.0_rl.fits')

hdulist_co12 = fits.open(cube_filename)
data_cube = hdulist_co12[0].data
wcs_co12 = WCS(hdulist_co12[0].header)

hdu_sig = fits.open(sig_file)
co12 = hdu_sig[0].data
wcs_co12 = WCS(hdu_sig[0].header)

def getPsi(path_to_file):
    I, Q, U, __, wcs = IQU(path_to_file)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, psi

I, __, __, wcs_250, psi = getPsi(blast250_file)

cdelta_blast = 0.00277777777778 * 60 # deg * 60 arcmin/deg
nskip_blast = np.int(np.round(1.0/cdelta_blast)) # points per arcmin
r_blast = nskip_blast*3.0

cdelta_co12 = 8.33333333333E-03 * 60 # deg * 60 arcmin/deg
nskip_co12 = np.int(np.round(1.0/cdelta_co12)) # points per arcmin
r_co12 = nskip_co12*3.0

# points to use for S calculation
xs = np.arange(psi.shape[1])[::nskip_blast]
xs = xs[30:-30] # don't use edges
ys = np.arange(psi.shape[0])[::nskip_blast]
ys = ys[5:-5] # don't use edges
X,Y = np.meshgrid(xs,ys)
mask_ra, mask_dec = wcs_250.all_pix2world(X, Y, 0)
mask_gal_coords = SkyCoord(mask_ra*u.degree, mask_dec*u.degree, frame='fk5').galactic

Ix = np.arange(I.shape[1])
Iy = np.arange(I.shape[0])
IX, IY = np.meshgrid(Ix, Iy)
ra, dec = wcs_250.all_pix2world(IX, IY, 0)
radec = SkyCoord(ra*u.degree, dec*u.degree, frame='fk5')

x_co12 = np.arange(co12.shape[1])
y_co12 = np.arange(co12.shape[0])
XCO12, YCO12 = np.meshgrid(x_co12, y_co12)
co12_ra, co12_dec = wcs_co12.all_pix2world(XCO12, YCO12, 0)
radec_co12 = SkyCoord(co12_ra*u.degree, co12_dec*u.degree, frame='galactic').fk5

#pol_disp = np.empty_like(I)
#pol_disp[:] = np.nan
pol_disp = np.load('pol_disp.npy')
pol_coords = np.column_stack((X.ravel(), Y.ravel()))
"""
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
    pol_disp[coord[1], coord[0]] = S
"""
co12_x, co12_y = wcs_co12.all_world2pix(mask_gal_coords.l, mask_gal_coords.b, 0)
# integer pixel values
co12_x = np.round(co12_x).astype('int')
co12_y = np.round(co12_y).astype('int')

co12_mask_ra, co12_mask_dec = wcs_co12.all_pix2world(co12_x, co12_y, 0)
co12_mask_fk5 = SkyCoord(co12_mask_ra*u.degree, co12_mask_dec*u.degree, frame='galactic').fk5

# make velocity dispersion map
#vdisp = np.empty_like(co12_x)
#vdisp[:] = np.nan

vcoords = np.column_stack((co12_x.ravel(), co12_y.ravel()))
"""
v_array = np.empty((len(vmask_coords), data_cube.shape[0]))
count = 0
for coord in vmask_coords:
    print count
    if coord[0] > 120 or coord[1] < 0:
        count += 1
        continue
    center = co12[coord[1], coord[0]]
    mask = (x_co12[np.newaxis,:] - coord[0])**2 + (y_co12[:,np.newaxis] - coord[1])**2 <= (r_co12)**2
    #mask3D = np.repeat(mask[np.newaxis,:, :,], data_cube.shape[0], axis=0)
    #masked_cube = np.ma.masked_array(data_cube, mask3D)
    #v_array[count] = masked_cube.mean(axis = (2,1))
    for i in range(data_cube.shape[0]):
        if not data_cube[i][mask].shape[0]:
            continue
        else:
            v_array[count][i] = np.nanmean(data_cube[i][mask])
            #print data_cube[i][mask].shape[0]
    count += 1
"""
vdisp = np.load('v_array.npy')
#vdisp_lpf, v_fwhm = fwhm(vdisp)
v_fwhm = np.load('fwhm.npy')
v_fwhm *= 1000.0
vdisp_map = np.empty_like(co12)
vdisp_map[:] = np.nan

count = 0
for vcoord in vcoords:
    if vcoord[0] > 120 or vcoord[1] < 0:
        count += 1
        continue
    vdisp_map[vcoord[1],vcoord[0]] = v_fwhm[count]
    count += 1

"""
B = np.empty_like(co12)
B[:] = np.nan
count = 0
for i in range(len(vcoords)):
    x = vcoords[i][0]
    y = vcoords[i][1]
    if x > 120 or y < 0:
        count += 1
        continue
    #B[y,x] = np.sqrt(I[pol_coords[i][1],pol_coords[i][0]]) * \
    #     (vdisp_map[y,x] / pol_disp[pol_coords[i][1],pol_coords[i][0]])
    #B[y,x] = I[pol_coords[i][1],pol_coords[i][0]]
    count += 1
"""

B = np.empty_like(pol_disp)
B[:] = np.nan
count = 0
for i in range(len(pol_coords)):
    y = pol_coords[i][1]
    x = pol_coords[i][0]
    vx = vcoords[i][0]
    vy = vcoords[i][1]
    if vx > 120 or vy < 0:
        count += 1
        continue
    B[y,x] = np.sqrt(I[y,x]) * (vdisp_map[vy,vx] / pol_disp[y,x])
    #B[y,x] = I[y,x]
    #B[y,x] = pol_disp[y,x]
    #B[y,x] = vdisp_map[vy,vx]

x, y = np.indices(B.shape)
B_interp = np.array(B)
B_interp[np.isnan(B_interp)] = griddata((x[~np.isnan(B_interp)], y[~np.isnan(B_interp)]),\
       B[~np.isnan(B_interp)],(x[np.isnan(B_interp)], y[np.isnan(B_interp)]))

plt.figure()
plt.xlim(np.max(np.asarray(radec.ra)), np.min(np.asarray(radec.ra)))
plt.ylim(np.min(np.asarray(radec.dec)), np.max(np.asarray(radec.dec)))
#plt.pcolormesh(radec.ra, radec.dec, I, cmap = "inferno")
#plt.pcolormesh(radec.ra, radec.dec, B, cmap = "inferno")
plt.contourf(radec.ra, radec.dec, I, levels = np.linspace(np.nanmin(I),\
          np.nanmax(I), 40), cmap = 'inferno', alpha = 0.8)
plt.contourf(radec.ra, radec.dec, B_interp,\
     levels = np.linspace(np.nanmin(B_interp),\
     np.nanmax(B_interp), 30), cmap = "viridis", alpha = 0.7)
plt.xlabel("RA (deg)")
plt.ylabel("DEC (deg)")
#plt.pcolormesh(radec.ra, radec.dec, pol_disp, alpha = 0.5)
#plt.contourf(radec.ra, radec.dec, B_interp,\
#     levels = np.linspace(np.nanmin(B_interp) + 500, np.nanmax(B_interp), 10), cmap = 'inferno')
#plt.pcolormesh(radec_co12.ra, radec_co12.dec, B, cmap = "inferno")
#plt.pcolormesh(radec_co12.ra, radec_co12.dec, vdisp_map, cmap = "viridis", alpha = 0.5)
#plt.pcolormesh(radec_co12.ra, radec_co12.dec, co12, cmap = "viridis", alpha = 0.5)
#plt.scatter(co12_mask_fk5.ra, co12_mask_fk5.dec, c = 'r', alpha = 0.5)
#plt.scatter(mask_ra, mask_dec, c = 'g', alpha = 0.5)

#plt.pcolormesh(XCO12, YCO12, co12, cmap = "viridis", alpha = 0.5)
#plt.scatter(co12_x, co12_y, c = 'r', alpha = 0.5)

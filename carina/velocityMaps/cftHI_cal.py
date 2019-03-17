import numpy as np
import matplotlib
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
from streamLines import plot_streams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import aplpy
plt.ion()
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"

info = np.loadtxt("./HI_INFO", dtype = "str")
root_dir = info[np.where(info == 'ROOT_DIR')[0][0]][1]
cube_file = os.path.join(root_dir,\
           info[np.where(info == 'CUBE_FILE')[0][0]][1])
blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST500_FILE')[0][0]][1])

vmin = float(info[np.where(info == 'VMIN')[0][0]][1])
vmax = float(info[np.where(info == 'VMAX')[0][0]][1])

cdelta_v = float(info[np.where(info == 'CDELTA_V')[0][0]][1])
cdelta_blast = float(info[np.where(info == 'CDELTA_BLAST')[0][0]][1])

#pol_disp = np.load('large_files/pol_disp_CO13_3.npy')
pol_disp = np.load('pol_disp_5arcmin.npy')
pol_disp = np.rad2deg(pol_disp)
v_array = np.load('./large_files/v_array_HI_2.npy')
v_fwhm = np.load('./large_files/fwhmHI_2.npy')
nH2 = np.load('nh2_est.npy')

hdulist_v = fits.open(cube_file)
data_cube = hdulist_v[0].data
dummy_v = hdulist_v[0].data[0]

wcs_v = WCS(naxis=2)
wcs_v.wcs.crpix = [563.0, 665.0]
wcs_v.wcs.cdelt = np.array([-0.00277777798786, 0.00277777798786])
wcs_v.wcs.crval = [287.41231676, -1.14784038739]
wcs_v.wcs.ctype = ["GLON-SIN", "GLAT-SIN"]

pol_eff = [0.81, 0.79, 0.82]

def getPsi(path_to_file, band_idx):
    I, Q, U, __, wcs = IQU(path_to_file)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    pvals[pvals > 0.5] = np.nan
    pvals[pvals < 0] = np.nan
    pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, pvals, psi

I, Q, U, wcs, pvals, psi = getPsi(blast_file, 2)

dx = pvals*np.cos(psi)
dy = pvals*np.sin(psi)

cdelta_blast *= 60 # deg * 60 arcmin/deg
nskip_blast = np.int(np.round(1.0/cdelta_blast)) # points per arcmin
r_blast = nskip_blast*1.5

cdelta_v *= 60 # deg * 60 arcmin/deg
nskip_v = np.int(np.round(1.0/cdelta_v)) # points per arcmin
r_v = nskip_v*1.5

# points to use for S calculation
xs = np.arange(psi.shape[1])[::nskip_blast]
xs = xs[30:-30] # don't use edges
ys = np.arange(psi.shape[0])[::nskip_blast]
ys = ys[5:-5] # don't use edges
X,Y = np.meshgrid(xs,ys)
vectors = np.array([dx,dy])

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

pol_coords = np.column_stack((X.ravel(), Y.ravel()))
"""
pol_disp = np.empty_like(I)
pol_disp[:] = np.nan
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
pol_disp = np.rad2deg(pol_disp)
"""
v_x, v_y = wcs_v.all_world2pix(mask_gal_coords.l, mask_gal_coords.b, 0)
# integer pixel values
v_x = np.round(v_x).astype('int')
v_y = np.round(v_y).astype('int')

v_mask_ra, v_mask_dec = wcs_v.all_pix2world(v_x, v_y, 0)
v_mask_fk5 = SkyCoord(v_mask_ra*u.degree, v_mask_dec*u.degree, frame='galactic').fk5

# make velocity dispersion map
vcoords = np.column_stack((v_x.ravel(), v_y.ravel()))
"""
v_array = np.empty((len(vcoords), data_cube.shape[0]))
count = 0
for coord in vcoords:
    print count
    if coord[0] > 120 or coord[1] < 0:
        count += 1
        continue
    center = dummy_v[coord[1], coord[0]]
    mask = (x_v[np.newaxis,:] - coord[0])**2 + (y_v[:,np.newaxis] - coord[1])**2 <= (r_v)**2
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

vdisp_lpf, v_fwhm = fwhm(v_array, vmin, vmax)
"""
vdisp_map = np.empty_like(dummy_v)
vdisp_map[:] = np.nan
count = 0
for vcoord in vcoords:
    if vcoord[0] > 120 or vcoord[1] < 0:
        count += 1
        continue
    vdisp_map[vcoord[1],vcoord[0]] = v_fwhm[count]
    count += 1

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
    # Calibrate into microgauss (uG) using Crutcher 2004 (Scuba) Eq 2
    # nH2 = g/cm^3, V = km/s
    B[y,x] = 9.3*np.sqrt(nH2[y,x]) * ((vdisp_map[vy,vx]/1.0e3) / pol_disp[y,x])
    #B[y,x] = (vdisp_map[vy,vx]/1.0e3) / pol_disp[y,x]

    #B[y,x] = np.sqrt(I[y,x]) * (vdisp_map[vy,vx] / pol_disp[y,x])
    #B[y,x] = I[y,x]
    #B[y,x] = pol_disp[y,x]
    #B[y,x] = vdisp_map[vy,vx]

x, y = np.indices(B.shape)
B_interp = np.array(B)
B_interp[np.isnan(B_interp)] = griddata((x[~np.isnan(B_interp)], y[~np.isnan(B_interp)]),\
       B[~np.isnan(B_interp)],(x[np.isnan(B_interp)], y[np.isnan(B_interp)]))

print "B max (G) = ", np.max(B_interp[np.isfinite(B_interp)])
print "B min (G) = ", np.min(B_interp[np.isfinite(B_interp)])
"""
plt.figure()
plt.xlim(np.max(np.asarray(radec.ra)), np.min(np.asarray(radec.ra)))
plt.ylim(np.min(np.asarray(radec.dec)), np.max(np.asarray(radec.dec)))
#plt.pcolormesh(radec.ra, radec.dec, I, cmap = "inferno")
#plt.pcolormesh(radec.ra, radec.dec, B, cmap = "inferno")
plt.contourf(radec.ra, radec.dec, I, levels = np.linspace(np.nanmin(I),\
          np.nanmax(I), 40), cmap = 'inferno', alpha = 0.8)
plt.contourf(radec.ra, radec.dec, B_interp,\
     levels = np.linspace(np.nanmin(B_interp),\
     np.nanmax(B_interp), 10), cmap = "viridis", alpha = 0.7)
plt.xlabel("RA (deg)")
plt.ylabel("DEC (deg)")
"""
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

hdu = fits.PrimaryHDU(data=B_interp, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('white')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)

dx = dx[30:-30,230:-270]
dy = dy[30:-30,230:-270]
IX = IX[30:-30,230:-270]
IY= IY[30:-30,230:-270]
vectors = np.array([dx,dy])
plot_streams(ax, vectors, IX, IY, nskip = 30, alph = 0.4, col = 'yellow', vec = False)

plt.imshow(B_interp, origin = 'lower', cmap = 'inferno', vmax = 300)
#im = ax.imshow(B, origin = 'lower', cmap = 'inferno', vmax = 50)

#plt.xlim(210, 660)
#plt.ylim(40, 440)

#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#cb1 = matplotlib.colorbar.ColorbarBase(cax, cmap='inferno', orientation='vertical').ax.tick_params(axis='y', direction='out')
#cax.xaxis.set_label_position('bottom')
#cax.xaxis.set_ticks_position('bottom')
#plt.colorbar(im, cax)

#plt.savefig(os.path.join(save_files_here, 'cftC013.png'), format='png', bbox_inches = 'tight')


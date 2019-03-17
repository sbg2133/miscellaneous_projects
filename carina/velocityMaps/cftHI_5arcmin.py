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
from vdispHI import fwhm as fwhm
from streamLines import plot_streams
from streamLines import plot_vectors
import os
import glob
from makePretty import pretty
from mpl_toolkits.axes_grid1 import make_axes_locatable
from fwhmHI import vfwhm
from columnDensity import nH2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ion()
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"
band = sys.argv[1] # '250', '350', or '500'
L = np.float(sys.argv[2]) # column depth, pc

info = np.loadtxt("./HI_INFO", dtype = "str")
root_dir = info[np.where(info == 'ROOT_DIR')[0][0]][1]
cube_file = os.path.join(root_dir,\
           info[np.where(info == 'CUBE_FILE')[0][0]][1])
if band == '250':
    blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST250_FILE')[0][0]][1])
if band == '350':
    blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST350_FILE')[0][0]][1])
if band == '500':
    blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST500_FILE')[0][0]][1])
cdelta_v = float(info[np.where(info == 'CDELTA_V')[0][0]][1])
cdelta_blast = float(info[np.where(info == 'CDELTA_BLAST')[0][0]][1])
Vmin = np.float(info[np.where(info == 'VMIN')[0][0]][1])
Vmax = np.float(info[np.where(info == 'VMAX')[0][0]][1])

res = 5.0 # arcmin

cdelta_blast *= 60 # arcmin * 60 arcmin/deg
nskip_blast = np.int(np.round(1.0/cdelta_blast)) # points per arcmin
r_blast = nskip_blast*(res/2.)

cdelta_v *= 60 # deg * 60 arcmin/deg
nskip_v = np.int(np.round(1.0/cdelta_v)) # points per arcmin
r_v = nskip_v*(res/2.)

pol_disp = np.load('pol_disp_skip3.npy')
v_array = np.load('./large_files/v_array_skip3.npy')
v_fwhm = np.load('./large_files/v_fwhm_nskip3.npy')
nskip_blast = 3

#pol_disp = np.load('pol_disp_default.npy')
#v_array = np.load('v_array_default.npy')

hdulist_v = fits.open(cube_file)
data_cube = hdulist_v[0].data
dummy_v = hdulist_v[0].data[0]

wcs_v = WCS(naxis=2)
wcs_v.wcs.crpix = [563.0, 665.0]
wcs_v.wcs.cdelt = np.array([-0.00277777798786, 0.00277777798786])
wcs_v.wcs.crval = [287.41231676, -1.14784038739]
wcs_v.wcs.ctype = ["GLON-SIN", "GLAT-SIN"]

pol_eff = [0.81, 0.79, 0.82]
cen_coord = [160.84429, -59.582361]

#I, Q, U, wcs, pvals, psi = getPsi(blast_file, 2)
I, Q, U, pol_data, wcs = IQU(blast_file, do_cov = True)
psi = np.deg2rad(pol_data[1])
sig_psi = np.deg2rad(pol_data[5])
pvals = pol_data[2]
pvals[pvals > 0.8] = np.nan
pvals[pvals < 0] = np.nan
dx = np.cos(psi)
dy = np.sin(psi)
n, N = nH2(band, I, L)

# points to use for S calculation
xs = np.arange(psi.shape[1])[::nskip_blast]
xs = xs[30:-30] # don't use edges
ys = np.arange(psi.shape[0])[::nskip_blast]
ys = ys[5:-5] # don't use edges
X,Y = np.meshgrid(xs,ys)
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
# Convert the HI galactic coordinates to fk5
radec_v = SkyCoord(v_ra*u.degree, v_dec*u.degree, frame='galactic').fk5

pol_coords = np.column_stack((X.ravel(), Y.ravel()))
sky_coords = np.column_stack((mask_ra.ravel(), mask_dec.ravel()))

"""
pol_disp = np.empty_like(I)
pol_disp[:] = np.nan

count = len(pol_coords)
for coord in pol_coords:
    print count
    center = psi[coord[1], coord[0]]
    var_center = sig_psi[coord[1],coord[0]]**2
    mask_inner = np.sqrt((Ix[np.newaxis,:] - coord[0])**2 + (Iy[:,np.newaxis] - coord[1])**2) < (r_blast)
    mask_outer = np.sqrt((Ix[np.newaxis,:] - coord[0])**2 + (Iy[:,np.newaxis] - coord[1])**2) <= (r_blast)
    mask = mask_inner ^ mask_outer
    comp_points = psi[mask]
    var_comp_points = sig_psi[mask]**2
    center = psi[coord[1], coord[0]]
    #if np.abs(np.rad2deg(center)) > 25.0:
    #    print "phi =", np.rad2deg(center)
    #    pol_disp[coord[1], coord[0]] = np.nan
    #    count += 1
    #    continue
    S2 = (center - comp_points)**2
    S2 = np.sum(S2)/len(comp_points)
    varS2 = np.sum(var_center + var_comp_points)/len(var_comp_points)
    S2_debias = S2 - varS2
    if S2_debias < 0:
        S2_debias = np.nan
    #S2_debias = np.sum(S2_db)/len(S2_db)
    #varS2 = np.sum(var_center + var_comp_points)
    #S2_debias = S2 - np.var(S2)
    #print S2, varS2, S2_debias
    #S = np.sqrt((1.0/len(S2))*np.sum(S2))
    #S = np.sqrt((1.0/len(S2_debias))*np.sum(S2_debias))
    S = np.sqrt(S2_debias)
    #I[coord[1], coord[0]] = 1.0e3
    pol_disp[coord[1], coord[0]] = S
    count -= 1
pol_disp = np.rad2deg(pol_disp)
"""
# Convert the HI galactic coordinates into pixel coords
v_x, v_y = wcs_v.all_world2pix(mask_gal_coords.l, mask_gal_coords.b, 0)
# integer pixel values
v_x = np.round(v_x).astype('int')
v_y = np.round(v_y).astype('int')

# Convert the pixel coordinates into fk5
v_mask_ra, v_mask_dec = wcs_v.all_pix2world(v_x, v_y, 0)
v_mask_fk5 = SkyCoord(v_mask_ra*u.degree, v_mask_dec*u.degree, frame='galactic').fk5

# make velocity dispersion map
vcoords = np.column_stack((v_x.ravel(), v_y.ravel()))
#v_array = np.empty((len(vcoords), data_cube.shape[0]))
steps = v_array.shape[1]
Vres = (Vmax - Vmin)/steps
V = np.arange(Vmin, Vmax, Vres)/1.0e3
"""
count = 0
for coord in vcoords:
    print len(vcoords) - count
    center = dummy_v[coord[1], coord[0]]
    mask = np.sqrt((x_v[np.newaxis,:] - coord[0])**2 + (y_v[:,np.newaxis] - coord[1])**2) <= (r_v)
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
print "Calculating v FWHM"
#v_fwhm = vfwhm(v_array)

def Bmap_Vmap(pol_coords, pol_disp, vcoords, v_fwhm):
    vdisp_map = np.empty_like(dummy_v)
    vdisp_map[:] = np.nan
    count = 0
    for vcoord in vcoords:
        vdisp_map[vcoord[1],vcoord[0]] = v_fwhm[count]
        count += 1
    B = np.empty_like(pol_disp)
    B[:] = np.nan
    V = np.empty_like(pol_disp)
    V[:] = np.nan
    v_sky_coords = np.column_stack((np.asarray(v_mask_fk5.ra).ravel(), np.asarray(v_mask_fk5.dec).ravel()))
    print "Calculating B"
    count = 0
    for i in range(len(pol_coords)):
        y = pol_coords[i][1]
        x = pol_coords[i][0]
        vy = vcoords[i][1]
        vx = vcoords[i][0]
        #print sky_coords[i][0], sky_coords[i][1]
        #print v_sky_coords[i][0], v_sky_coords[i][1]
        #print vy, vx
        #print
        # Calibrate into microgauss (uG) using Crutcher 2004 (Scuba) Eq 2
        # n = g/cm^3, V = km/s
        B[y,x] = 0.2*9.3*np.sqrt(2.0*n[y,x]) * (vdisp_map[vy,vx] /\
                  (np.sqrt(2.0*pol_disp[y,x]**2)) )
        #B[y,x] = np.sqrt(I[y,x]) * (vdisp_map[vy,vx] / pol_disp[y,x])
        #B[y,x] = I[y,x]
        #B[y,x] = pol_disp[y,x]
        V[y,x] = vdisp_map[vy,vx]
    x, y = np.indices(B.shape)
    V_interp = np.array(V)
    V_interp[np.isnan(V_interp)] = griddata((x[~np.isnan(V_interp)], y[~np.isnan(V_interp)]),\
           V[~np.isnan(V_interp)],(x[np.isnan(V_interp)], y[np.isnan(V_interp)]))
    print "V max (km/s) = ", np.max(V_interp[np.isfinite(V_interp)])
    print "V min (km/s) = ", np.min(V_interp[np.isfinite(V_interp)])
    print
    B_interp = np.array(B)
    B_interp[np.isnan(B_interp)] = griddata((x[~np.isnan(B_interp)], y[~np.isnan(B_interp)]),\
           B[~np.isnan(B_interp)],(x[np.isnan(B_interp)], y[np.isnan(B_interp)]))
    print "B max (uG) = ", np.max(B_interp[np.isfinite(B_interp)])
    print "B min (uG) = ", np.min(B_interp[np.isfinite(B_interp)])
    return B, B_interp, V, V_interp

B, B_interp, V, V_interp = Bmap_Vmap(pol_coords, pol_disp, vcoords, v_fwhm)

mask_coords = np.column_stack((v_x.ravel(), v_y.ravel()))
for coord in mask_coords[::120]:
#    if coord[0] > 120 or coord[1] < 0:
#        continue
    center = dummy_v[coord[1], coord[0]]
    mask = np.sqrt((x_v[np.newaxis,:] - coord[0])**2 + (y_v[:,np.newaxis] - coord[1])**2) <= (r_v)
    comp_points = dummy_v[mask]
    dummy_v[mask] = 1.0e8

"""
plt.figure(figsize = (12,12))
disp = np.mean(vdisp, axis = 0)
V = V[:-1]
plt.plot(V, disp, linewidth = 2)
ax = plt.gca()
ax.set_xlabel('Radial Velocity [km/s]')
ax.set_ylabel('Counts')
pretty()
plt.savefig(os.path.join(save_files_here, 'HI_avg_vdisp.eps'), format='eps', bbox_inches = 'tight')
"""
###########################################################
# HISTOGRAM of V FWHM
###########################################################
plt.figure(figsize = (12,12))
plt.hist(v_fwhm, bins = 500, histtype='bar', ec='black', color = 'C0')
plt.xlabel(r"FWHM V$_{LOS}$ [km/s]")
plt.ylabel('N')
plt.xlim(14, 39)
pretty()
plt.savefig(os.path.join(save_files_here, 'VLOS_hist.eps'), format='eps', bbox_inches = 'tight')

##########################################################
# HISTOGRAM OF NH2 Column density
###########################################################
plt.figure(figsize = (12,12))
plt.hist(np.log10(N), bins = 500, histtype='bar', ec='black', color = 'C0')
plt.xlabel(r"log$_{10}$ N(H$_{2})$ cm$^{-2}$")
plt.ylabel('N')
plt.xlim(19.8, 22.8)
pretty()
plt.savefig(os.path.join(save_files_here, 'NH2_hist_' + str(band) + '.eps'),\
          format='eps', bbox_inches = 'tight')

#########################################################
# V select regions
#########################################################
plt.figure(figsize = (12,12))
plt.pcolormesh(radec_v.ra, radec_v.dec, dummy_v, cmap = "viridis", alpha = 0.7)
plt.contourf(radec.ra, radec.dec, I, levels = np.linspace(np.nanmin(I),\
          np.nanmax(I), 40), cmap = 'inferno', alpha = 0.5)
ax = plt.gca()
ax.set_ylim(-58.5, -60.5)
ax.set_xlim(158, 164)
ax.set_ylabel('RA (J2000)')
ax.set_xlabel('Dec (J2000)')
ax.set_xlim(ax.get_xlim()[::-1])
pretty()
plt.savefig(os.path.join(save_files_here, 'HI_vselect.png'),\
          format='png', bbox_inches = 'tight')


###########################################################
# B interp with pol streamlines
###########################################################
hdu = fits.PrimaryHDU(data=B_interp, header=wcs.to_header())
f = aplpy.FITSFigure(hdu, figsize = (12, 12))
f.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno', smooth = 3,\
     interpolation = 'hanning')
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
f.ticks.set_color('white')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.show()
f.colorbar.set_axis_label_text(r'$B_{pos}$ $\mu$G')
f.colorbar.set_axis_label_font(size = 16)

dx = dx[30:-30,230:-270]
dy = dy[30:-30,230:-270]
IX = IX[30:-30,230:-270]
IY= IY[30:-30,230:-270]
vectors = np.array([dx,dy])
#plot_vectors(ax, vectors, IY, IX, nskip = 30, alph = 0.4, col = 'white', pot = False)
plot_streams(ax, vectors, IX, IY, nskip = 30, alph = 0.4, col = 'yellow', vec = False)
plt.savefig(os.path.join(save_files_here, 'cftHI_5arcmin_' + str(band)\
             + '_' + str(L) + 'pc' + '.png'), format='png', bbox_inches = 'tight')

"""
hdu = fits.PrimaryHDU(data=V_interp, header=wcs.to_header())
fv = aplpy.FITSFigure(hdu, figsize = (12, 12))
fv.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
fv.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
fv.show_colorscale(cmap = 'viridis', smooth = 3,\
     interpolation = 'hanning')
fv.add_scalebar(15/60.) # arcmin
fv.scalebar.set_color('white')
fv.scalebar.set_corner('bottom right')
fv.scalebar.set_label('10 pc')
fv.scalebar.set_linewidth('2')
fv.scalebar.set_font_size('16')
fv.tick_labels.set_yformat('dd.dd')
fv.tick_labels.set_xformat('dd.dd')
fv.axis_labels.set_font(size=16)
fv.tick_labels.set_font(size=16)
fv.ticks.set_color('white')
fv.ticks.set_linewidth('2')
fv.add_colorbar()
fv.colorbar.set_location('right')
fv.colorbar.show()
fv.colorbar.set_axis_label_text(r'$\sigma_{v}$ [km/s]')
fv.colorbar.set_axis_label_font(size = 16)
"""
"""
# H2 Number density (1/cm^3)
hdu = fits.PrimaryHDU(data=nH2, header=wcs.to_header())
f = aplpy.FITSFigure(hdu, figsize = (12, 12))
f.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno', smooth = 3,\
     interpolation = 'hanning')
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
f.ticks.set_color('white')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.show()
f.colorbar.set_axis_label_text(r'cm$^{-3}$')

dx = dx[30:-30,230:-270]
dy = dy[30:-30,230:-270]
IX = IX[30:-30,230:-270]
IY= IY[30:-30,230:-270]
vectors = np.array([dx,dy])
#plot_vectors(ax, vectors, IY, IX, nskip = 30, alph = 0.4, col = 'white', pot = False)
plot_streams(ax, vectors, IX, IY, nskip = 30, alph = 0.5, col = 'white', vec = False)
#plt.savefig(os.path.join(save_files_here, 'nH2_5arcmin.png'), format='png', bbox_inches = 'tight')

dx = dx[30:-30,230:-270]
dy = dy[30:-30,230:-270]
IX = IX[30:-30,230:-270]
IY= IY[30:-30,230:-270]
"""
V = V[30:-30,230:-270]
B = B[30:-30,230:-270]
I = I[30:-30,230:-270]
#vectors = np.array([dx,dy])

###############################
# B LIC FIGURE
###############################
lic = np.loadtxt("../lic.dat")
lic = np.transpose(lic)
mult = lic*B_interp[30:-30,260:-260]
mult -= np.nanmin(mult)
plt.figure(figsize = (12,12))
plt.subplot(projection=wcs)
ax = plt.gca()
plt.imshow(mult, cmap = 'inferno', vmin = 70, vmax  = 110, interpolation = 'bilinear')
#ax.set_facecolor("k")
plt.xlabel('RA (J2000)')
plt.ylabel('Dec (J2000)')
#ax.coords['ra'].set_axislabel('RA (J2000)')
#ax.coords['dec'].set_axislabel('Dec (J2000)')
pretty()
plt.savefig(os.path.join(save_files_here, 'B_lic_' + str(band)\
       + '_' + str(L) + 'pc' + '.eps'), format='eps', bbox_inches = 'tight')

################################################
# V FWHM scatter with pol vectors
################################################
plt.figure(figsize = (12,12))
plt.subplot(projection=wcs)
plt.scatter(IX, IY, c = V, cmap = 'viridis', vmin = 15, vmax = 35)
ax = plt.gca()
#ax.set_facecolor("k")
plot_vectors(ax, vectors, IY, IX, nskip = 20, alph = 0.3, col = 'white', pot = False)
ax.coords['ra'].set_axislabel('RA (J2000)')
ax.coords['dec'].set_axislabel('Dec (J2000)')
ax.tick_params(axis='both', width = 2)
cbar = plt.colorbar(pad = 0.01)
cbar.set_label(r'V$_{FWHM}$ [km/s]')
plt.savefig(os.path.join(save_files_here, 'VFWHM_vectors.png'), format='png', bbox_inches = 'tight')

####################################
# B scatter with pol vectors
###################################
plt.figure(figsize = (12,12))
plt.subplot(projection=wcs)
#ax.set_facecolor("k")
plt.scatter(IX, IY, c = B, cmap = 'inferno', vmin = 10, vmax = 400)
ax = plt.gca()
cbar = plt.colorbar(pad = 0.01)
cbar.set_label(r'B$_{pos}$ [$\mu$G]')
ax.coords['ra'].set_axislabel('RA (J2000)')
ax.coords['dec'].set_axislabel('Dec (J2000)')
ax.tick_params(axis='both', width = 2)
plot_vectors(ax, vectors, IY, IX, nskip = 20, alph = 0.3, col = 'white', pot = False)
plt.savefig(os.path.join(save_files_here, 'B_vectors_' + str(band)\
              + '_' + str(L) + 'pc' + '.png'), format='png', bbox_inches = 'tight')

#overlay = ax.get_coords_overlay('galactic')
#overlay[0].set_axislabel('Galactic Longitude')
#overlay[1].set_axislabel('Galactic Latitude')
#overlay.grid(color='red', linestyle='solid', alpha=0.5)
#ax.coords['ra'].set_ticks(color='white')
#ax.coords['dec'].set_ticks(color='white')
#ax.coords['ra'].set_axislabel('RA (J2000)')
#ax.coords['dec'].set_axislabel('Dec (J2000)')
#ax.coords.grid(color='red', linestyle='solid', alpha=0.5)

B0p5pc = np.load('Bhist_500_0.5pc.npy')
B5pc = np.load('Bhist_500_5pc.npy')
B5pc = B5pc[B5pc > 15.0]
B0p5pc = B0p5pc[B0p5pc > 40.0]
plt.figure(figsize = (12,12))
plt.hist(B0p5pc, bins = 500, histtype='bar', ec='black',\
      color = 'C0', alpha = 0.7, label = 'L = 0.5 pc')
plt.hist(B5pc, bins = 500, histtype='bar', ec='black',\
      color = 'purple', alpha = 0.7, label = 'L = 5 pc')
plt.xlabel(r"B$_{pos}$ $\mu$G")
plt.ylabel('N')
plt.xlim(10, 900)
plt.savefig(os.path.join(save_files_here, 'B_hist_' + str(band) + '.eps'),\
       format='eps', bbox_inches = 'tight')
plt.legend(loc = 'upper right')
pretty()

n5pc = np.load('nH2_500_5pc_hist.npy')
n0p5pc = np.load('nH2_500_0p5pc_hist.npy')
fig = plt.figure(figsize = (17,12))
ax = fig.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax1 = fig.add_subplot(121)
plt.hist(n0p5pc, bins = 500, histtype='bar', ec='black',\
            color = 'C0', alpha = 0.7, label = 'L = 0.5 pc')
ax1.set_xlim(0, 6000)
plt.legend(loc = 'upper right')
pretty()
ax2 = fig.add_subplot(122, sharey = ax1)
ax2.axes.get_yaxis().set_visible(False)
plt.hist(n5pc, bins = 500, histtype='bar', ec='black',\
            color = 'purple', alpha = 0.7, label = 'L = 5 pc')
ax2.set_xlim(0, 600)
plt.legend(loc = 'upper right')
pretty()
ax.set_xlabel(r"n(H$_{2}$) cm$^{-3}$")
ax.set_ylabel('N')
fig.subplots_adjust(wspace=0,hspace=0)
plt.savefig(os.path.join(save_files_here, 'nH2_hist_' + str(band) + '.png'),\
       format='png', bbox_inches = 'tight')

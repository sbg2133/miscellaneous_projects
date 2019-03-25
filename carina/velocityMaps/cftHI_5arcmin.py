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
#from vdispHI import fwhm as fwhm
from streamLines import plot_streams
from streamLines import plot_vectors
import os
import glob
from makePretty import pretty
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from fwhmHI import vfwhm
from fitFWHM import vfwhm
from columnDensity import nH2
from add_regions import car_regions as reg
from scipy import stats
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ion()
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"

# plot?
plot = 1

# save plots?
save = 1

# debias p?
debias = 1

# calc S?
calc_S = 0

band = sys.argv[1] # '250', '350', or '500'
L = np.float(sys.argv[2]) # column depth, pc
l = str(L)
if '.' in l:
    l = l.replace('.', 'p')
fit = sys.argv[3] # 'lor', 'gauss'

info = np.loadtxt("./HI_INFO", dtype = "str")
root_dir = info[np.where(info == 'ROOT_DIR')[0][0]][1]
cube_file = os.path.join(root_dir,\
           info[np.where(info == 'CUBE_FILE')[0][0]][1])
if band == '250':
    blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST250_FILE')[0][0]][1])
    band_idx = 0
if band == '350':
    blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST350_FILE')[0][0]][1])
    band_idx = 1
if band == '500':
    blast_file = os.path.join(root_dir,\
           info[np.where(info == 'BLAST500_FILE')[0][0]][1])
    band_idx = 2
cdelta_v = float(info[np.where(info == 'CDELTA_V')[0][0]][1])
cdelta_blast = float(info[np.where(info == 'CDELTA_BLAST')[0][0]][1])
Vmin = np.float(info[np.where(info == 'VMIN')[0][0]][1])
Vmax = np.float(info[np.where(info == 'VMAX')[0][0]][1])

res = 5.0 # arcmin
mH2 = 3.32e-24

cdelta_blast *= 60 # arcmin * 60 arcmin/deg
nskip_blast = np.int(np.round(1.0/cdelta_blast)) # points per arcmin
r_blast = nskip_blast*(res/2.)

cdelta_v *= 60 # deg * 60 arcmin/deg
nskip_v = np.int(np.round(1.0/cdelta_v)) # points per arcmin
r_v = nskip_v*(res/2.)

if not calc_S:
    pol_disp = np.load('pol_disp_nskip3_3.npy')
v_array = np.load('./large_files/v_array_skip3.npy')
#v_fwhm = np.load('./large_files/v_fwhm_nskip3.npy')
if fit == 'lor':
    v_fwhm = np.load('./large_files/fwhm_HI_lor.npy')
else:
    v_fwhm = np.load('./large_files/fwhm_HI_gauss.npy')
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
if debias:
    I, Q, U, pol_data, wcs = IQU(blast_file, do_cov = True)
    psi = pol_data[1] # deg
    sig_psi = pol_data[4] # deg
    sig_p = pol_data[5]
    pvals = pol_data[2]
    #pvals[pvals > 0.5] = np.nan
    #pvals[pvals < 0] = np.nan
    dx = np.cos(np.deg2rad(psi))
    dy = np.sin(np.deg2rad(psi))
    # Galactic plane = 27 deg
    #dx = np.cos(np.deg2rad((27.0)*np.ones_like(psi)))
    #dy = np.sin(np.deg2rad((27.0)*np.ones_like(psi)))
    #dx = np.cos(np.deg2rad(27.0 + psi))
    #dy = np.sin(np.deg2rad(27.0 + psi))
else:
    I, Q, U, __, wcs = IQU(blast_file, do_cov = False)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    pvals[pvals > 0.5] = np.nan
    pvals[pvals < 0] = np.nan
    pvals /= pol_eff[band_idx]
    psi = np.rad2deg(0.5*np.arctan2(U,Q))
    dx = np.cos(np.deg2rad(psi))
    dy = np.sin(deg2rad(psi))

# H2 Number density and column density
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

if calc_S:
    pol_disp = np.empty_like(I)
    pol_disp[:] = np.nan
    count = len(pol_coords)
    for coord in pol_coords:
        print count
        center = psi[coord[1], coord[0]]
        mask_inner = np.sqrt((Ix[np.newaxis,:] - coord[0])**2 + (Iy[:,np.newaxis] - coord[1])**2) < (r_blast)
        mask_outer = np.sqrt((Ix[np.newaxis,:] - coord[0])**2 + (Iy[:,np.newaxis] - coord[1])**2) <= (r_blast)
        mask = mask_inner ^ mask_outer
        comp_points = psi[mask]
        center = psi[coord[1], coord[0]]
        #if np.abs(np.rad2deg(center)) > 25.0:
        #    print "phi =", np.rad2deg(center)
        #    pol_disp[coord[1], coord[0]] = np.nan
        #    count += 1
        #    continue
        S = center - comp_points
        S2 = np.nansum(S**2)/len(S[~np.isnan(S)])
        S2 -= np.nanvar(S2)
        Sdb = np.sqrt(S2)
        print Sdb
        #if debias:
        #    varS2 = np.nanvar(S)
            #var_center = sig_psi[coord[1],coord[0]]**2
            #var_comp_points = sig_psi[mask]**2
            #varS2 = np.nansum(var_center + var_comp_points)/len(var_comp_points[~np.isnan(var_comp_points)])
        #else:
        #    varS2 = np.nanvar(S)
        #S2_debias = S2 - varS2
        #Sdb = np.sqrt(S2_debias)
        #print S2, varS2, Sdb

        #I[coord[1], coord[0]] = 1.0e3
        pol_disp[coord[1], coord[0]] = Sdb
        count -= 1

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
#vel, v_fwhm = vfwhm(v_array)
print "<v fwhm> =", np.nanmean(v_fwhm)
print "std v dist =", np.nanstd(v_fwhm)

want_coords = []
center = [160.69,-59.59]
want_coords.append(center)
sb_r = [160.37,-59.78]
want_coords.append(sb_r)
sb_l = [160.88,-59.75]
want_coords.append(sb_l)
sb_bright = [160.61,-60.08]
want_coords.append(sb_bright)
tr_16 = [161.41,-59.74]
want_coords.append(tr_16)
tr_14 = [160.95,-59.55]
want_coords.append(tr_14)
nb = [160.56,-59.42]
want_coords.append(nb)
sp = [161.76,-60.01]
want_coords.append(sp)
want_coords = np.asarray(want_coords)

tab_coords = []
coord_idxs = []
for coord in want_coords:
    dist = np.sum((sky_coords - coord)**2, axis=1)
    tab_coords.append(sky_coords[np.argmin(dist)])
    coord_idxs.append(np.argmin(dist))
tab_coords = np.asarray(tab_coords)
# Vals are: ra, dec, Bpos, nH2, S, vfwhm
tab_vals = np.zeros((len(tab_coords),6))

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
    S = np.empty_like(pol_disp)
    S[:] = np.nan
    #v_sky_coords = np.column_stack((np.asarray(v_mask_fk5.ra).ravel(), np.asarray(v_mask_fk5.dec).ravel()))
    print "Calculating B"
    count = 0
    for i in range(len(pol_coords)):
        # add values to table array for tab_coords
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
        V[y,x] = vdisp_map[vy,vx]
        S[y,x] = pol_disp[y,x]
        if i in coord_idxs:
            #print i, sky_coords[i][0]
            tab_idx = np.where(sky_coords[i][0] == tab_coords[:,0])[0][0]
            # print tab_idx
            # ra
            tab_vals[tab_idx][0] = tab_coords[tab_idx][0]
            # dec
            tab_vals[tab_idx][1] = tab_coords[tab_idx][1]
            # Bpos
            tab_vals[tab_idx][2] = B[y,x]
            # nH2
            tab_vals[tab_idx][3] = n[y,x]
            # S
            tab_vals[tab_idx][4] = S[y,x]
            # Vfwhm
            tab_vals[tab_idx][5] = V[y,x]

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
    return B, B_interp, V, V_interp, S

B, B_interp, V, V_interp, S = Bmap_Vmap(pol_coords, pol_disp, vcoords, v_fwhm)

b = B[~np.isnan(B)]
print "<B> =", np.mean(b)
print "sigma B =", np.std(b)
np.save('tab_vals.npy', tab_vals)

############################################################################
# Southern bubble B estimate (see: https://iopscience.iop.org/article/10.1088/0004-637X/760/2/150/pdf
###########################################################################

# add_regions.py
Rout = 43 # arcmin
#Rout = 43 # arcmin
rout = 41 # arcmin
Rin = 28 # arcmin
# rin = 20 # arcmin
rin = 26 # arcmin

# Churchwell, Eq. 2
R = 0.5*(np.sqrt(Rout*rout) + np.sqrt(Rin*rin))
T = np.sqrt(Rout*rout) - np.sqrt(Rin*rin)

T_frac = T/R
eta = 0.29
cs = 10.0 # km/s (T = 10^4)
cs = 1.0e4 * 1.0e2 # cm / s

amu_g = 1.66e-24 # g
#rho_ism = 1.4 * amu_g # g cm^-3, diffuse ISM
rho_ism = 10.0 * mH2 # g

# Pavel 2012, Eq. 5

# rho_shell / rho_ism = 1 / (1 - T)**3
rho_rat = 1.0/ (1.0 - ((1.0 - T_frac)**3.0))

rho_shell = rho_rat * rho_ism
#rho_shell = 1.78 * rho_ism

B_est1 = 1.0e6 * np.sqrt(0.5 * rho_shell * cs**2.0 * 8.0 * np.pi)
B_est2 = 1.0e6 * np.sqrt(0.5 * eta * rho_shell * (cs**2.0) * 8.0 * np.pi)
print "B est1 lower limit (uG) =", B_est1
print "B est2 lower limit (uG) =", B_est2

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
if save:
    plt.savefig(os.path.join(save_files_here, 'HI_avg_vdisp.eps'), format='eps', bbox_inches = 'tight')
"""
"""
dx = dx[30:-30,230:-270]
dy = dy[30:-30,230:-270]
IX = IX[30:-30,230:-270]
IY= IY[30:-30,230:-270]
vectors = np.array([dx,dy])

plt.figure()
ax = plt.subplot(projection = wcs)
plot_vectors(ax, vectors, IY, IX, nskip = 20, alph = 0.4, col = 'k', pot = False)
overlay = ax.get_coords_overlay('galactic')
overlay.grid(color='red', linestyle='solid', alpha=0.5)
"""
Phi = psi[30:-30,260:-260]
z = Phi + 127.0 # 90 deg is parallel to galactic plane
z = z[~np.isnan(z)]
z[z > 180.0] -= 90
print "<z> (deg) =", np.mean(z)
print "sigma z (deg) =", np.std(z)
# Percentage of vectors within +- 15 deg of galactic plane
perc = len(z[np.where(np.logical_and(z>=(90 - 23.), z<=(90.0 + 23)))[0]])/float(len(z))
print "% within +- 23 deg of gal plane =", perc
Smean = np.mean(S[~np.isnan(S)])
print "<S> (deg) =", np.mean(S[~np.isnan(S)])
print "sigma S (deg) =", np.std(S[~np.isnan(S)])


if plot:
    ##########################################################
    # Angle between Phi and galactic plane (90 deg = parallel)
    plt.figure(figsize = (12,12))
    plt.hist(z, bins = 50, histtype='bar', ec='black', color = 'C0')
    plt.xlabel(r"$\Phi$ (deg)")
    plt.ylabel('Number of Sightlines')
    pretty()
    if save:
        plt.savefig(os.path.join(save_files_here, 'Phi_hist.eps'), format='eps', bbox_inches = 'tight')

    ###########################################################
    # HISTOGRAM of S
    ###########################################################
    S = S[30:-30,230:-270]
    S = S[~np.isnan(S)]
    S = S[S > 0]
    plt.figure(figsize = (12,12))
    plt.hist(S, bins = 50, histtype='bar', ec='black', color = 'C0')
    plt.xlabel(r"S$_{\Phi}$ [deg]")
    plt.ylabel('Number of Sightlines')
    plt.xlim(0, 100)
    pretty()
    if save:
        plt.savefig(os.path.join(save_files_here, 'S_hist.eps'), format='eps', bbox_inches = 'tight')

    ###########################################################
    # HISTOGRAM of V FWHM
    ###########################################################
    Vtemp = V[30:-30,230:-270]
    Vtemp = Vtemp[~np.isnan(Vtemp)]
    Vtemp = Vtemp[Vtemp > 0]
    plt.figure(figsize = (12,12))
    plt.hist(Vtemp, bins = 50, histtype='bar', ec='black', color = 'C0')
    plt.xlabel(r"FWHM V$_{LOS}$ [km/s]")
    plt.ylabel('Number of Sightlines')
    if fit == 'lor':
        plt.xlim(20, 57)
    else:
        plt.xlim(23, 50)
    pretty()
    if save:
        plt.savefig(os.path.join(save_files_here, 'VLOS_hist.eps'), format='eps', bbox_inches = 'tight')

    ##########################################################
    # HISTOGRAM OF NH2 Column density
    ###########################################################
    plt.figure(figsize = (12,12))
    plt.hist(np.log10(N), bins = 100, histtype='bar', ec='black', color = 'C0')
    plt.xlabel(r"log$_{10}$ N(H$_{2})$ cm$^{-2}$")
    plt.ylabel('Number of Sightlines')
    plt.xlim(19.8, 22.8)
    pretty()
    if save:
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
    if save:
        plt.savefig(os.path.join(save_files_here, 'HI_vselect.png'),\
              format='png', bbox_inches = 'tight')

    ###########################################################
    # B interp with pol streamlines
    ###########################################################
    hdu = fits.PrimaryHDU(data=B_interp, header=wcs.to_header())
    f = aplpy.FITSFigure(hdu)
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
    plot_streams(ax, vectors, IX, IY, nskip = 30, alph = 0.5, col = 'yellow', vec = False)
    reg(ax, f, wcs, c = 'lime', ls = '--', add_dots = True)
    if save:
        plt.savefig(os.path.join(save_files_here, 'cftHI_5arcmin_' + str(band)\
                 + '_' + l + 'pc' + '.png'), format='png', bbox_inches = 'tight')

    V = V[30:-30,230:-270]
    B = B[30:-30,230:-270]
    I = I[30:-30,230:-270]

    ###############################
    # B LIC FIGURE
    ###############################
    wcs.wcs.crpix = np.array([-40.0, 25.0])
    lic = np.loadtxt("../lic.dat")
    lic -= np.nanmin(lic)
    lic = np.transpose(lic)
    mult = lic*B_interp[30:-30,260:-260]
    mult -= np.nanmin(mult)
    hdu = fits.PrimaryHDU(data=mult, header=wcs.to_header())
    f = aplpy.FITSFigure(hdu)
    f.set_theme('publication')
    ax = plt.gca()
    ax.set_facecolor("k")
    f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
    f.show_colorscale(cmap = 'inferno', vmin = 1, vmax = 400, interpolation = 'hanning')
    #f.add_scalebar(15/60.) # arcmin
    #f.scalebar.set_color('white')
    #f.scalebar.set_corner('bottom right')
    #f.scalebar.set_label('10 pc')
    #f.scalebar.set_linewidth('2')
    #f.scalebar.set_font_size('16')
    f.tick_labels.set_yformat('dd.dd')
    f.tick_labels.set_xformat('dd.dd')
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size=16)
    f.ticks.set_color('white')
    f.ticks.set_linewidth('2')
    reg(ax, f, wcs, c = 'lime', ls = '--', add_dots = True)
    wcs.wcs.crpix = np.array([ 220.,   50.])
    #f.add_colorbar()
    #f.colorbar.set_location('right')
    #f.colorbar.show()
    #f.colorbar.set_axis_label_text(r'$B_{pos}$ $\mu$G')
    #f.colorbar.set_axis_label_font(size = 16)
    """
    plt.figure(figsize = (12,12))
    plt.subplot(projection=wcs)
    ax = plt.gca()
    if fit == 'lor':
        plt.imshow(mult, cmap = 'inferno', origin = 'lower', vmin = 1, vmax  = 350, interpolation = 'bilinear')
    else:
        plt.imshow(mult, cmap = 'inferno', origin = 'lower', vmin = 1, vmax  = 400, interpolation = 'bilinear')
    #ax.set_facecolor("k")
    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)')
    #ax.coords['ra'].set_axislabel('RA (J2000)')
    #ax.coords['dec'].set_axislabel('Dec (J2000)')
    pretty()
    """
    if save:
        plt.savefig(os.path.join(save_files_here, 'B_lic_' + str(band)\
           + '_' + l + 'pc' + '.eps'), format='eps', bbox_inches = 'tight')

    ################################################
    # V FWHM scatter with pol vectors
    ################################################
    plt.figure(figsize = (12,12))
    plt.subplot(projection=wcs)
    if fit == 'lor':
        plt.scatter(IX, IY, c = V, cmap = 'viridis', vmin = 20)
    else:
        plt.scatter(IX, IY, c = V, cmap = 'viridis', vmin = 30)
    ax = plt.gca()
    #ax.set_facecolor("k")
    plot_vectors(ax, vectors, IY, IX, nskip = 20, alph = 0.3, col = 'white', pot = False)
    ax.coords['ra'].set_axislabel('RA (J2000)')
    ax.coords['dec'].set_axislabel('Dec (J2000)')
    ax.tick_params(axis='both', width = 2)
    cbar = plt.colorbar(pad = 0.01)
    cbar.set_label(r'V$_{FWHM}$ [km/s]')
    if save:
        plt.savefig(os.path.join(save_files_here, 'VFWHM_vectors.png'), format='png', bbox_inches = 'tight')

    def colorbar(mappable):
        ax = mappable.axes
        fig = ax.figure
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        return fig.colorbar(mappable, cax=cax)

    ####################################
    # B scatter with pol vectors
    ###################################
    fig = plt.figure(figsize = (12, 14))
    plt.subplot(projection=wcs)
    ax = plt.gca()
    ax.set_facecolor("k")
    if fit == 'lor':
        plt.scatter(IX, IY, c = B, cmap = 'inferno', vmin = 10, vmax = 720)
    else:
        plt.scatter(IX, IY, c = B, cmap = 'inferno', vmin = 10, vmax = 720)
    #cbar = plt.colorbar(pad = 0.1)
    #cbar = fig.colorbar(b, orientation="horizontal", pad=0.1, aspect = 40)
    #cbar.set_label(r'B$_{pos}$ [$\mu$G]')
    ax.coords['ra'].set_axislabel('RA (J2000)')
    ax.coords['dec'].set_axislabel('Dec (J2000)')
    #ax.tick_params(axis='both', which = 'both',  width = 2, color = 'white')
    plot_vectors(ax, vectors, IY, IX, nskip = 20, alph = 0.3, col = 'white', pot = False)
    overlay = ax.get_coords_overlay('galactic')
    overlay.grid(color='red', linestyle='solid', alpha=0.8)
    reg(ax, f, wcs, c = 'lime', ls = '--', add_dots = True)
    if save:
        plt.savefig(os.path.join(save_files_here, 'B_vectors_' + str(band)\
                  + '_' + l + 'pc' + '.png'), format='png', bbox_inches = 'tight')

    #ax.coords['ra'].set_ticks(color='white')
    #ax.coords['dec'].set_ticks(color='white')
    #ax.coords['ra'].set_axislabel('RA (J2000)')
    #ax.coords['dec'].set_axislabel('Dec (J2000)')
    #ax.coords.grid(color='red', linestyle='solid', alpha=0.5)

    #############################################################
    # B HISTOGRAM for 5 pc and 0.5 pc
    #################################################################
    B0p5pc = np.load('Bhist_500_0p5pc.npy')
    B5pc = np.load('Bhist_500_5pc.npy')
    B5pc = B5pc[B5pc > 0.0]
    B0p5pc = B0p5pc[B0p5pc > 0.0]
    plt.figure(figsize = (12,12))
    #plt.hist(B0p5pc, bins = 150, histtype='bar', ec='black',\
    #      color = 'C0', alpha = 0.7, label = 'L = 0.5 pc')
    plt.hist(B5pc, bins = 100, histtype='bar', ec='black',\
          color = 'purple', alpha = 0.7, label = 'L = 5 pc')
    plt.xlabel(r"B$_{pos}$ $\mu$G")
    plt.ylabel('Number of Sightlines')
    plt.xlim(0, 560)
    plt.legend(loc = 'upper right')
    pretty()
    if save:
        plt.savefig(os.path.join(save_files_here, 'B_hist_' + str(band) + '.eps'),\
           format='eps', bbox_inches = 'tight')
    ##################################################
    # nH2 HISTOGRAM for 5 pc
    ###################################################
    n5pc = np.load('nH2_500_5pc_hist.npy')
    n0p5pc = np.load('nH2_500_0p5pc_hist.npy')
    """
    fig = plt.figure(figsize = (17,12))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax1 = fig.add_subplot(121)
    plt.hist(n0p5pc, bins = 100, histtype='bar', ec='black',\
                color = 'C0', alpha = 0.7, label = 'L = 0.5 pc')
    ax1.set_xlim(0, 6000)
    plt.legend(loc = 'upper right')
    pretty()
    ax2 = fig.add_subplot(122, sharey = ax1)
    ax2.axes.get_yaxis().set_visible(False)
    plt.hist(n5pc, bins = 150, histtype='bar', ec='black',\
                color = 'purple', alpha = 0.7, label = 'L = 5 pc')
    ax2.set_xlim(0, 600)
    plt.legend(loc = 'upper right')
    pretty()
    ax.set_xlabel(r"n(H$_{2}$) cm$^{-3}$")
    ax.set_ylabel('Number of Sightlines')
    fig.subplots_adjust(wspace=0,hspace=0)
    """
    plt.figure(figsize = (12,12))
    plt.hist(n5pc, bins = 200, histtype='bar', ec='black',\
          color = 'purple', alpha = 0.7, label = 'L = 5 pc')
    plt.xlabel(r"n(H$_{2}$) cm$^{-3}$")
    plt.ylabel('Number of Sightlines')
    plt.xlim(0, 600)
    plt.legend(loc = 'upper right')
    pretty()
    if save:
        plt.savefig(os.path.join(save_files_here, 'nH2_hist_' + str(band) + '.eps'),\
           format='eps', bbox_inches = 'tight')


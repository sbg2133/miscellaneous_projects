import numpy as np
import matplotlib.pyplot as plt
from getIQU import IQU
from scipy import ndimage
from astropy import units as u
from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from streamLines import plot_vectors
from streamLines import plot_streams
import aplpy

##################################################################
# run like: plotIntensity(['250'], streamlines = True, nskip = 12)
##################################################################
plt.ion()
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']
pol_eff = [0.81, 0.79, 0.82]
all_bands = ['250', '350', '500']
cen_coord = [160.84429, -59.582361]

im250 = './carinaData/carinaneb/carinaneb_good_250_p10_good_C_gls_map_cal.fits'
im350 = './carinaData/carinaneb/carinaneb_good_350_p10_good_C_gls_map_cal.fits'
im500 = './carinaData/carinaneb/carinaneb_good_500_p10_good_C_gls_map_cal.fits'
smooth250 = './carinaData/smooth/3.0_arcmin/carinaneb_250_smoothed_3.0_rl.fits'
smooth350 = './carinaData/smooth/3.0_arcmin/carinaneb_350_smoothed_3.0_rl.fits'
smooth500 = './carinaData/smooth/3.0_arcmin/carinaneb_500_smoothed_3.0_rl.fits'
planck850 = './carinaData/planckData/planck_353_carinaneb_pol.fits'
smooth500_4p8um = './carinaData/smooth/4.8_arcmin/carinaneb_500_smoothed_4.8_rl.fits'

#I250, Q250, U250, wcs = IQU('250', im250)
#I350, Q350, U350, wcs = IQU('350', im350)
#I500, Q500, U500, wcs = IQU('500', im500)

#I250, Q250, U250, wcs = IQU('250', smooth250)
#I350, Q350, U350, wcs = IQU('350', smooth350)
#I500, Q500, U500, wcs = IQU('500', smooth500)

#Is = [I250, I350, I500]
#Qs = [Q250, Q350, Q500]
#Us = [U250, U350, U500]

"""
I250[np.isnan(I350)] = 0.
I250[I250 < 0.] = 0.
I350[np.isnan(I350)] = 0.
I350[I350 < 0.] = 0.
I500[np.isnan(I500)] = 0.
I500[I500 < 0.] = 0.
I250[I250 < 5*np.std(I250)] = 0.
I350[I350 < 5*np.std(I350)] = 0.
I500[I500 < 5*np.std(I500)] = 0.
"""

# FK5 Eta Carina coords (SIMBAD): 10 45 03.546 -59 41 03.95
EC_dec = Angle('-59d41m3.95s')
EC_ra = Angle('10h45m3.546s')
# Keyhole Nebula 10 44 19.0 -59 53 21
KH_dec = Angle('-59d53m21s')
KH_ra = Angle('10h44m19.0s')
# WR25 10 44 10.391 -59 43 11.10
WR_dec = Angle('-59d43m11.10s')
WR_ra = Angle('10h44m10.391s')
# HD93205 10 44 33.739 -59 44 15.44
HD93_dec = Angle('-59d44m15.44s')
HD93_ra = Angle('10h44m33.739s')
# Trumpler 16 10 45 10.0 -59 43 0
T16_dec = Angle('-59d43m0s')
T16_ra = Angle('10h45m10s')
# Trumpler 14 10 43 56  -59 33 00
TR14_dec = Angle('-59d33m0s')
TR14_ra = Angle('10h43m56s')
# CPD 59.2661, star in "Treasure Chest" 10 46 54 -59.2661
CPD_dec = Angle('-59d57m0s')
CPD_ra = Angle('10h46m54s')
marker_dec = np.array([EC_dec.deg, WR_dec.deg, HD93_dec.deg, CPD_dec.deg])
marker_ra = np.array([EC_ra.deg, WR_ra.deg, HD93_ra.deg, CPD_ra.deg])
labels = ['eta Car', 'WR25', 'HD93205', 'CPD'] 

def getPsi(filename):
    I, Q, U, wcs = IQU(filename)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, psi

def addScalebar(fig, label):
    #  scalebar
    fig.add_scalebar(15/60.) # arcmin
    fig.scalebar.set_color('white')
    fig.scalebar.set_corner('bottom right')
    fig.scalebar.set_label(label)
    fig.scalebar.set_linewidth(2)
    fig.scalebar.set_font_size(14)
    return

def plotIntensity(filename, streamlines = False, vec = True, nskip = 20, annotate = False):
    #plt.style.use('dark_background')
    f = aplpy.FITSFigure(filename)
    f.set_theme('publication')
    plt.tight_layout()
    f.show_colorscale(cmap = 'inferno')
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size = 14)
    addScalebar(f, '10 pc')
    f.add_grid()
    f.grid.set_color('yellow')
    f.grid.set_alpha(0.3)
    ra, dec = marker_ra, marker_dec
    if annotate:
        f.show_markers(ra, dec, edgecolor='white', facecolor='none',
                marker='o', s=100, alpha=0.5)
        for i in range(len(labels)):
            f.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
    if (vec) or (streamlines):
        ax = plt.gca()
        I, Q, U, __, psi = getPsi(filename)
        dx = np.cos(psi)
        dy = np.sin(psi)
        X = np.linspace(0, I.shape[1], I.shape[1])
        Y = np.linspace(0, I.shape[0], I.shape[0])
        xs, ys = np.meshgrid(X,Y)
        vectors = np.array([dx,dy])
    if vec:
        plot_vectors(ax, vectors, ys, xs, nskip = nskip, alph = 0.4, col = 'white')
    if streamlines:
        plot_streams(ax, vectors, xs, ys, nskip = nskip, alph = 0.2, col = 'yellow', vec = False)
        #putStreamlines(ax, smooth250, nskip, vec = True, alph = 0.5, col = 'yellow')
        #putStreamlines(ax, planck850, nskip, vec = False, alph = 0.5, col = 'orange')
        #putStreamlines(ax, smooth350, nskip, vec = True, alph = 0.5, col = 'orange')
        #putStreamlines(ax, smooth500, nskip, vec = True, alph = 0.5, col = 'red')
    ax.set_facecolor("k")
    plt.tight_layout()
    plt.savefig('./intensity.png', dpi = 100, bbox_inches = 'tight')
    return

def overplotPlanck(filename, streamlines = False, vec = True, nskip = 20, annotate = False):
    fnames = [filename, planck850]
    vec_colors = ['white', 'yellow']
    f = aplpy.FITSFigure(filename)
    #f.set_theme('publication')
    plt.tight_layout()
    f.show_colorscale(cmap = 'inferno')
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size = 14)
    addScalebar(f, '10 pc')
    f.add_grid()
    f.grid.set_color('yellow')
    f.grid.set_alpha(0.3)
    ra, dec = marker_ra, marker_dec
    if annotate:
        f.show_markers(ra, dec, edgecolor='white', facecolor='none',
                marker='o', s=100, alpha=0.5)
        for i in range(len(labels)):
            f.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
    for i in range(len(fnames)):
        if (vec) or (streamlines):
            ax = plt.gca()
            I, Q, U, __, psi = getPsi(fnames[i])
            dx = np.cos(psi)
            dy = np.sin(psi)
            X = np.linspace(0, I.shape[1], I.shape[1])
            Y = np.linspace(0, I.shape[0], I.shape[0])
            xs, ys = np.meshgrid(X,Y)
            vectors = np.array([dx,dy])
        if vec:
            plot_vectors(ax, vectors, ys, xs, nskip = nskip, alph = 0.4, col = vec_colors[i])
        if streamlines:
            plot_streams(ax, vectors, xs, ys, nskip = nskip, alph = 0.2, col = vec_colors[i], vec = False)
    ax.set_facecolor("k")
    plt.tight_layout()
    plt.savefig('./planck_500um_vec.png', dpi = 100, bbox_inches = 'tight')
    return

def overplotBands(filename, streamlines = False, vec = True, nskip = 20, annotate = False):
    fnames = [filename, smooth350, smooth500]
    vec_colors = ['white', 'yellow', 'red']
    f = aplpy.FITSFigure(filename)
    #f.set_theme('publication')
    plt.tight_layout()
    f.show_colorscale(cmap = 'inferno')
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size = 14)
    addScalebar(f, '10 pc')
    f.add_grid()
    f.grid.set_color('yellow')
    f.grid.set_alpha(0.3)
    ra, dec = marker_ra, marker_dec
    if annotate:
        f.show_markers(ra, dec, edgecolor='white', facecolor='none',
                marker='o', s=100, alpha=0.5)
        for i in range(len(labels)):
            f.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
    for i in range(len(fnames)):
        if (vec) or (streamlines):
            ax = plt.gca()
            I, Q, U, __, psi = getPsi(fnames[i])
            dx = np.cos(psi)
            dy = np.sin(psi)
            X = np.linspace(0, I.shape[1], I.shape[1])
            Y = np.linspace(0, I.shape[0], I.shape[0])
            xs, ys = np.meshgrid(X,Y)
            vectors = np.array([dx,dy])
        if vec:
            plot_vectors(ax, vectors, ys, xs, nskip = nskip, alph = 0.5, col = vec_colors[i])
        if streamlines:
            plot_streams(ax, vectors, xs, ys, nskip = nskip, alph = 0.2, col = vec_colors[i], vec = False)
    ax.set_facecolor("k")
    plt.tight_layout()
    plt.savefig('./overplot_bands_vec.png', dpi = 100, bbox_inches = 'tight')
    return

def gradDots(filename, annotate = False, streamlines = False, nskip = 20, plot = False, zoom = False):
    I, __, __, wcs, psi = getPsi(filename)
    I = I[30:-30,260:-260]
    #I += np.abs(np.nanmin(I))
    psi = psi[30:-30,260:-260]
    gradi, gradq = ndimage.sobel(I), ndimage.sobel(I, axis = 0)
    i = np.cos(psi)
    q = np.sin(psi)
    gradi, gradq, i, q, = np.nan_to_num(gradi), np.nan_to_num(gradq), np.nan_to_num(i), np.nan_to_num(q)
    ii = gradi * i
    iq = gradi * q
    qi = gradq * i
    qq = gradq * q
    print np.mean(ii + qq)/np.std(ii + qq)
    print np.mean(iq + qi)/np.std(iq + qi)
    #print 'sum gradi*i =', ii
    #print 'sum gradi*q =', iq
    #print 'sum gradq*q =', qq
    #print 'sum gradq*i =', qi
    #plt.imshow(gradI*i, cmap = 'inferno', origin = 'lower')
    return gradi, gradq, i, q

def getCurl(filename, annotate = False, streamlines = False, nskip = 20, plot = False, zoom = False):
    __, __, __, wcs, psi = getPsi(filename)
    x = np.cos(psi)
    y = np.sin(psi)
    #curl = np.diff(y)[:-1,:] - np.diff(x, axis = 0)[:,:-1]
    curl = ndimage.sobel(psi, axis = 0) - ndimage.sobel(psi, axis = 1)
    if plot:
        hdu1 = fits.PrimaryHDU(data=curl, header=wcs.to_header())
        f1 = aplpy.FITSFigure(hdu1)
        f1.set_theme('publication')
        ax = plt.gca()
        ax.set_facecolor("k")
        f1.show_colorscale(cmap = 'inferno')
        f1.axis_labels.set_font(size=16)
        f1.tick_labels.set_font(size = 14)
        #  scalebar
        f1.add_scalebar(15/60.) # arcmin
        f1.scalebar.set_label('10 pc')
        f1.scalebar.set_color('white')
        f1.scalebar.set_corner('bottom right')
        f1.scalebar.set_label('10 pc')
        f1.scalebar.set_linewidth(2)
        f1.scalebar.set_font_size(size = 'large')
        f1.add_grid()
        f1.grid.set_color('yellow')
        f1.grid.set_alpha(0.3)
        f1.add_colorbar()
        f1.colorbar.set_location('right')
        f1.colorbar.show()
        if zoom:
            f1.recenter(cen_coord[0], cen_coord[1], width = 1.5, height = 1.25)
        ra, dec = marker_ra, marker_dec
        if annotate:
            f1.show_markers(ra, dec, edgecolor='white', facecolor='none',
                    marker='o', s=100, alpha=0.5)
            for i in range(len(labels)):
                f1.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
        if streamlines:
            ax = plt.gca()
            putStreamlines(ax, filename, nskip = 20)
        plt.tight_layout()
        plt.savefig('./curl.png', dpi = 100, bbox_inches = 'tight')
    return curl, wcs

def curlIntensity(filename, annotate = False, streamlines = False, nskip = 20, zoom = False):
    I, __, __, wcs, __ = getPsi(filename)
    curl, wcs = getCurl(smooth250)
    hdu1 = fits.PrimaryHDU(np.sqrt(curl**2), header=wcs.to_header())
    fig = plt.figure()
    f1 = aplpy.FITSFigure(hdu1, figure = fig)
    f1.set_theme('publication')
    ax = plt.gca()
    ax.set_facecolor("k")
    f1.show_colorscale(cmap = 'inferno')
    f1.axis_labels.set_font(size=16)
    f1.tick_labels.set_font(size = 14)
    #  scalebar
    f1.add_scalebar(15/60.) # arcmin
    f1.scalebar.set_label('10 pc')
    f1.scalebar.set_color('white')
    f1.scalebar.set_corner('bottom right')
    f1.scalebar.set_label('10 pc')
    f1.scalebar.set_linewidth(2)
    f1.scalebar.set_font_size(size = 'large')
    f1.add_grid()
    f1.grid.set_color('yellow')
    f1.grid.set_alpha(0.3)
    f1.add_colorbar()
    f1.colorbar.set_location('right')
    f1.colorbar.show()
    if zoom:
        f1.recenter(cen_coord[0], cen_coord[1], width = 1.5, height = 1.25)
    ra, dec = marker_ra, marker_dec
    if annotate:
        f1.show_markers(ra, dec, edgecolor='white', facecolor='none',
                marker='o', s=100, alpha=0.5)
        for i in range(len(labels)):
            f1.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
    if streamlines:
        #ax = plt.gca()
        putStreamlines(ax, filename, nskip = 30, col = 'white', alph = 0.3)
    plt.tight_layout()
    #hdu1.writeto('./mag_curl.fits')
    plt.savefig('./curl_v_intensity.png', dpi = 100, bbox_inches = 'tight')
    return

def plotSobelCurl(filename, annotate = False, streamlines = False, nskip = 20, zoom = False):
    curl, wcs = getCurl(smooth250)
    #sx, sy = ndimage.sobel(curl, 0), ndimage.sobel(curl, axis = 1)
    #mag = np.sqrt(sx**2 + sy**2)
    hdu1 = fits.PrimaryHDU(data=np.sqrt(curl**2), header=wcs.to_header())
    fig = plt.figure()
    f1 = aplpy.FITSFigure(hdu1, figure = fig)
    f1.set_theme('publication')
    ax = plt.gca()
    ax.set_facecolor("k")
    f1.show_colorscale(cmap = 'inferno')
    f1.axis_labels.set_font(size=16)
    f1.tick_labels.set_font(size = 14)
    #  scalebar
    f1.add_scalebar(15/60.) # arcmin
    f1.scalebar.set_label('10 pc')
    f1.scalebar.set_color('white')
    f1.scalebar.set_corner('bottom right')
    f1.scalebar.set_label('10 pc')
    f1.scalebar.set_linewidth(2)
    f1.scalebar.set_font_size(size = 'large')
    f1.add_grid()
    f1.grid.set_color('yellow')
    f1.grid.set_alpha(0.3)
    f1.add_colorbar()
    f1.colorbar.set_location('right')
    f1.colorbar.show()
    if zoom:
        f1.recenter(cen_coord[0], cen_coord[1], width = 1.5, height = 1.25)
    ra, dec = marker_ra, marker_dec
    if annotate:
        f1.show_markers(ra, dec, edgecolor='white', facecolor='none',
                marker='o', s=100, alpha=0.5)
        for i in range(len(labels)):
            f1.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
    if streamlines:
        putStreamlines(ax, filename, nskip = 30, col = 'white', alph = 0.3, vec = True)
    plt.tight_layout()
    #hdu1.writeto('./sobel_curl.fits')
    plt.savefig('./sobel_curl.png', dpi = 100, bbox_inches = 'tight')
    return

def putStreamlines(ax, filename, nskip, alph = 1, col = 'yellow', vec = False):
    I, Q, U, __, psi = getPsi(filename)
    dx = np.cos(psi)
    dy = np.sin(psi)
    mag = np.sqrt(dx**2 + dy**2)
    X = np.linspace(0, I.shape[1], I.shape[1])
    Y = np.linspace(0, I.shape[0], I.shape[0])
    xs, ys = np.meshgrid(X,Y)
    xsize, ysize = len(X), len(Y)
    vectors = np.array([dx,dy])
    plot_streams(ax, vectors, xs, ys, nskip, col, alph, vec)
    return

#overplotIntensity(smooth500_4p8um, streamlines = False, vec = True, nskip = 30, annotate = True)
#plotIntensity(smooth250, streamlines = True, vec = True, nskip = 30, annotate = True)
#overplotBands(smooth250, streamlines = False, vec = True, nskip = 30, annotate = True)

#plotIntensity(smooth350, streamlines = True, nskip = 30, annotate = True)
#plotIntensity(smooth500, streamlines = True, nskip = 30, annotate = True)
#getCurl(smooth250, annotate = True, streamlines = False, plot = True, zoom = True)
#curlIntensity(smooth250, annotate = True, streamlines = True, zoom = False)
#plotSobelCurl(smooth250, annotate = True, streamlines = True, zoom = False)
#gradi, gradq, i, q = gradDots(smooth250, annotate = True, streamlines = True, zoom = False)

import numpy as np
import matplotlib.pyplot as plt
from getIQU import IQU
from streamLines import plot_streams
from scipy import ndimage
from astropy import units as u
from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
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
# CPD 59.2661, star in "Treasure Chest" 10 46 54 -59.2661
CPD_dec = Angle('-59d57m0s')
CPD_ra = Angle('10h46m54s')
marker_dec = np.array([EC_dec.deg, KH_dec.deg, WR_dec.deg, HD93_dec.deg, CPD_dec.deg])
marker_ra = np.array([EC_ra.deg, KH_ra.deg, WR_ra.deg, HD93_ra.deg, CPD_ra.deg])
labels = ['eta Car', 'Keyhole', 'WR25', 'HD93205', 'CPD'] 

def compositeImage():
    # Read in the three images downloaded from here:
    g = fits.open(im350)[1].data
    r = fits.open(im500)[1].data
    i = fits.open(im250)[1].data
    #Is[0][np.isnan(Is[0])] = 0.
    #Is[1][np.isnan(Is[1])] = 0.
    #Is[2][np.isnan(Is[2])] = 0.
    #Is[1] /= np.max(Is[1])
    #Is[2] /= np.max(Is[2])
    #Is[0] /= np.max(Is[0])
    rgb = make_lupton_rgb(i, r, g, Q=30, stretch=1000)
    plt.imshow(rgb, origin='lower')
    #img = aplpy.FITSFigure('cube_2d.fits')
    # Let's add a scalebar to it
    #img.add_scalebar(5/60.)
    #img.scalebar.set_label('5 arcmin')
    #img.scalebar.set_color('white')
    # We may want to lengthen the scalebar, move it to the top left,
    # and apply a physical scale
    #img.scalebar.set_corner('top left')
    #img.scalebar.set_length(17/60.)
    #img.scalebar.set_label('1 parsec')
    return

def getPsi(band_idx, filename):
    I, Q, U, wcs = IQU(band_idx, filename)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, psi

def plotIntensity(band_idx, filename, streamlines = False, nskip = 10, annotate = False, bdir = False):
    band_idx = 0
    #plt.style.use('dark_background')
    f = aplpy.FITSFigure(filename)
    f.set_theme('publication')
    plt.tight_layout()
    f.show_colorscale(cmap = 'inferno')
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size = 14)
    #  scalebar
    f.add_scalebar(30/60.) # arcmin
    f.scalebar.set_label('0.5 deg')
    f.scalebar.set_color('white')
    f.scalebar.set_corner('bottom right')
    f.scalebar.set_label('0.5 deg')
    f.scalebar.set_linewidth(2)
    f.scalebar.set_font_size(14)
    f.add_grid()
    f.grid.set_color('yellow')
    f.grid.set_alpha(0.3)
    ra, dec = marker_ra, marker_dec
    if annotate:
        f.show_markers(ra, dec, edgecolor='white', facecolor='none',
                marker='o', s=100, alpha=0.5)
        for i in range(len(labels)):
            f.add_label(ra[i] - 0.01, dec[i] - 0.01, labels[i], color='white')
    if streamlines:
        ax = plt.gca()
        putStreamlines(ax, filename, band_idx, nskip, bdir)
    plt.tight_layout()
    plt.savefig('./intensity.png', dpi = 100, bbox_inches = 'tight')
    return

def getCurl(band_idx, filename, annotate = False, streamlines = False, nskip = 20, bdir = False, plot = False, zoom = False):
    __, __, __, wcs, psi = getPsi(band_idx, filename)
    if bdir:
        psi += np.pi/4.
    x = np.cos(psi)
    y = np.sin(psi)
    curl = np.diff(y)[:-1,:] - np.diff(x, axis = 0)[:,:-1]
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
        f1.add_scalebar(30/60.) # arcmin
        f1.scalebar.set_label('0.5 deg')
        f1.scalebar.set_color('white')
        f1.scalebar.set_corner('bottom right')
        f1.scalebar.set_label('0.5 deg')
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
            putStreamlines(ax, filename, band_idx, nskip = 20)
        plt.tight_layout()
        plt.savefig('./curl.png', dpi = 100, bbox_inches = 'tight')
    return curl, wcs

def plotSobelCurl(band_idx, filename, annotate = False, streamlines = False, nskip = 20, bdir = False, zoom = False):
    curl, wcs = getCurl(band_idx, smooth250)
    sx, sy = ndimage.sobel(curl), ndimage.sobel(curl, axis = 0)
    hdu1 = fits.PrimaryHDU(data=np.sqrt(sx**2 + sy**2), header=wcs.to_header())
    f1 = aplpy.FITSFigure(hdu1)
    f1.set_theme('publication')
    ax = plt.gca()
    ax.set_facecolor("k")
    f1.show_colorscale(cmap = 'inferno')
    f1.axis_labels.set_font(size=16)
    f1.tick_labels.set_font(size = 14)
    #  scalebar
    f1.add_scalebar(30/60.) # arcmin
    f1.scalebar.set_label('0.5 deg')
    f1.scalebar.set_color('white')
    f1.scalebar.set_corner('bottom right')
    f1.scalebar.set_label('0.5 deg')
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
        putStreamlines(ax, filename, band_idx, nskip = 20)
    plt.tight_layout()
    plt.savefig('./sobel_curl.png', dpi = 100, bbox_inches = 'tight')
    return

def putStreamlines(ax, filename, band_idx, nskip, bdir = False):
    I, Q, U, __, psi = getPsi(band_idx, filename)
    if bdir:
        psi += np.pi/4.
    dx = np.cos(psi)
    dy = np.sin(psi)
    mag = np.sqrt(dx**2 + dy**2)
    X = np.linspace(0, I.shape[1], I.shape[1])
    Y = np.linspace(0, I.shape[0], I.shape[0])
    xs, ys = np.meshgrid(X,Y)
    xsize, ysize = len(X), len(Y)
    vectors = np.array([dx,dy])
    plot_streams(ax, vectors, xs, ys, nskip, alph = 0.2, col = 'yellow')
    return

#plot_streams(vectors, xs, ys, nskip = 10, vec = True)
#plotIntensity(0, smooth250, streamlines = True, nskip = 20, annotate = True)
getCurl(0, smooth250, annotate = True, bdir = True, streamlines = False, plot = True, zoom = True)
plotSobelCurl(0, smooth250, annotate = True, bdir = True, streamlines = False, zoom = True)

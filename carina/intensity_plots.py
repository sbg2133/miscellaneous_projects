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

im250 = './carinaData/carinaneb/carinaneb_good_250_p10_good_C_gls_map_cal.fits'
im350 = './carinaData/carinaneb/carinaneb_good_350_p10_good_C_gls_map_cal.fits'
im500 = './carinaData/carinaneb/carinaneb_good_500_p10_good_C_gls_map_cal.fits'
smooth250 = './carinaData/smooth/2.5_arcmin/carinaneb_250_smoothed_2.5_rl.fits'
smooth350 = './carinaData/smooth/2.5_arcmin/carinaneb_350_smoothed_2.5_rl.fits'
smooth500 = './carinaData/smooth/2.5_arcmin/carinaneb_500_smoothed_2.5_rl.fits'
#I250, Q250, U250, wcs = IQU('250', im250)
#I350, Q350, U350, wcs = IQU('350', im350)
#I500, Q500, U500, wcs = IQU('500', im500)

I250, Q250, U250, wcs = IQU('250', smooth250)
I350, Q350, U350, wcs = IQU('350', smooth350)
I500, Q500, U500, wcs = IQU('500', smooth500)

Is = [I250, I350, I500]
Qs = [Q250, Q350, Q500]
Us = [U250, U350, U500]

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

def plotIntensity(band, streamlines = False, nskip = 10, annotate = False):
    I = Is[bands.index(band)]
    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
    fig.add_axes(ax)
    ax.imshow(I, cmap = 'inferno', origin = 'lower')
    ax.coords.grid(color='yellow', alpha = 0.5, ls='solid')
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_axislabel('RA', fontsize = 14)
    dec.set_axislabel('DEC', fontsize = 14)
    dec.set_major_formatter('dd:mm:ss.s')
    ra.set_major_formatter('hh:mm:ss.s')
    ax.scatter(CPD_ra.deg, CPD_dec.deg, transform=ax.get_transform('world'), s=100, edgecolor='yellow', facecolor='none')
    ax.scatter(EC_ra.deg, EC_dec.deg, transform=ax.get_transform('world'), s=100, edgecolor='yellow', facecolor='none')
    ax.scatter(WR_ra.deg, WR_dec.deg, transform=ax.get_transform('world'), s=100, edgecolor='yellow', facecolor='none')
    ax.scatter(HD93_ra.deg, HD93_dec.deg, transform=ax.get_transform('world'), s=100, edgecolor='yellow', facecolor='none')
    ax.set_xlim(-0.5, I.shape[1] - 0.5)
    ax.set_ylim(-0.5, I.shape[0] - 0.5)
    overlay = ax.get_coords_overlay('galactic')
    overlay.grid(color='white', ls='dotted')
    overlay[0].set_axislabel('Galactic Longitude')
    overlay[1].set_axislabel('Galactic Latitude')
    if annotate:
        # eta carinae
        pix = wcs.wcs_world2pix(EC_ra.deg, EC_dec.deg, 0)
        ax.text(EC_ra.deg + 0.02, EC_dec.deg + 0.02,r'$\eta$Car', horizontalalignment='center',\
                         verticalalignment='center', color = 'white',\
                         fontsize = 10, transform=ax.get_transform('world'))
        # keyhole nebula
        ax.text(KH_ra.deg, KH_dec.deg, 'Keyhole', horizontalalignment='center',\
                         verticalalignment='center', color = 'white',\
                         fontsize = 10, transform=ax.get_transform('world'))
        # WR25
        ax.text(WR_ra.deg + 0.02, WR_dec.deg + 0.02, 'WR25', horizontalalignment='center',\
                         verticalalignment='center', color = 'white',\
                         fontsize = 10, transform=ax.get_transform('world'))
        # HD93
        ax.text(HD93_ra.deg + 0.02, HD93_dec.deg + 0.02, 'HD93205', horizontalalignment='center',\
                         verticalalignment='center', color = 'white',\
                         fontsize = 10, transform=ax.get_transform('world'))
        # T16
        ax.text(T16_ra.deg + 0.02, T16_dec.deg + 0.02, 'T16', horizontalalignment='center',\
                         verticalalignment='center', color = 'white',\
                         fontsize = 10, transform=ax.get_transform('world'))
        # CPD 59.2661
        ax.text(CPD_ra.deg + 0.02, CPD_dec.deg + 0.02, 'CPD 59.2661', horizontalalignment='center',\
                         verticalalignment='center', color = 'white',\
                         fontsize = 10, transform=ax.get_transform('world'))
    if streamlines:
        putStreamlines(ax, bands.index(band), nskip)
    return

def plotIntensityGrad(bands_list):
    f, ax = plt.subplots(len(bands_list), 1, sharex = True, squeeze=False, subplot_kw = {'projection': wcs})
    f.text(0.5, 0.04, 'RA', ha='center', fontsize = 12)
    f.text(0.04, 0.5, 'DEC', va='center', rotation='vertical', fontsize = 12)
    [ax[bands.index(band), 0].imshow(ndimage.sobel(Is[bands.index(band)]), cmap = "inferno") for band in bands_list]
    plt.tight_layout()
    return

def plotIntensityPsi(band):
    band_idx = bands.index(band)
    #I = Is[band_idx]
    #Q = Qs[band_idx]
    #U = Us[band_idx]
    I = Is[band_idx][:,260:-260]
    Q = Qs[band_idx][:,260:-260]
    U = Us[band_idx][:,260:-260]
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    #x = np.cos(psi)
    #y = np.sin(psi)
    #field = x + y
    sobel_curl = ndimage.sobel(psi, 1) - ndimage.sobel(psi, 0)
    f, ax = plt.subplots(3, sharex = True, subplot_kw = {'projection': wcs})
    curl = np.diff(psi, 1) - np.diff(psi, 0)[:,:-1]
    ax[0].imshow(curl, cmap = 'inferno')
    ax[1].imshow(ndimage.sobel(curl), cmap = 'inferno')
    ax[2].imshow(sobel_curl, cmap = 'inferno')
    #ax.imshow(ndimage.sobel(psi), cmap = 'inferno')
    #f.text(0.5, 0.04, 'RA', ha='center', fontsize = 12)
    #f.text(0.04, 0.5, 'DEC', va='center', rotation='vertical', fontsize = 12)
    #plt.tight_layout()
    return

def putStreamlines(ax, band_idx, nskip):
    """
    I = Is[band_idx][30:-30,260:-260]
    Q = Qs[band_idx][30:-30,260:-260]
    U = Us[band_idx][30:-30,260:-260]
    """
    I = Is[band_idx]
    Q = Qs[band_idx]
    U = Us[band_idx]
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # Correct pvals as in Jamil's thesis, 5.7
    #pvals[pvals > 0.5] = np.nan
    pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    dx = np.cos(psi)
    dy = np.sin(psi)
    mag = np.sqrt(dx**2 + dy**2)
    X = np.linspace(0, I.shape[1], I.shape[1])
    Y = np.linspace(0, I.shape[0], I.shape[0])
    xs, ys = np.meshgrid(X,Y)
    xsize, ysize = len(X), len(Y)
    vectors = np.array([dx,dy])
    plot_streams(ax, vectors, xs, ys, nskip)
    return

#plot_streams(vectors, xs, ys, nskip = 10, vec = True)
plotIntensity('250', streamlines = True, nskip = 20, annotate = True)


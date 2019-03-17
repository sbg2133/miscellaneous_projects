import os, sys
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from getIQU import IQU
from scipy import ndimage
from astropy import units as u
from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import aplpy
import matplotlib.image as mpimg
from streamLines import plot_vectors
from streamLines import plot_streams
#from intensity_plots_aplpy import getPsi
plt.ion()

save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"

# carina distance
d = 2.3e3
scale = d * 0.5 * np.pi/180.
smooth250 = './carinaData/smooth/3.0_arcmin_cal/carinaneb_250_smoothed_3.0_rl.fits'
smooth500 = './carinaData/smooth/3.0_arcmin_cal/carinaneb_500_smoothed_3.0_rl.fits'

overplot_dir = './carinaData/overplots'

#msx_png = os.path.join(overplot_dir, 'msx_250.png')
msx_png = os.path.join(overplot_dir, 'msx2.png')
hst_png = os.path.join(overplot_dir, 'hstgsc2_250.png')
irac_png = os.path.join(overplot_dir, 'irac_250.png')
spires_png = os.path.join(overplot_dir, 'spires250.png')
planck_png = os.path.join(overplot_dir, 'planck500.png')
curlIntensity_png = os.path.join(overplot_dir, 'sobel_curl2.png')
chandra_soft_png = os.path.join(overplot_dir, 'chandra_soft.png')
shassa32_png = os.path.join(overplot_dir, 'shassa32.png')

#files = [msx_png, hst_png, irac_png, spires_png, planck_png, curlIntensity_png]
files = [msx_png, shassa32_png, chandra_soft_png, hst_png, spires_png]

# load in I, Q, U for desired band
Ivals, Qvals, Uvals, __, wcs = IQU(smooth250)

I = Ivals
Q = Qvals
U = Uvals
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

# Correct pvals as in Jamil's thesis, 5.7
pvals[pvals > 0.5] = np.nan
pvals[pvals < 0] = np.nan
pol_eff = [0.81, 0.79, 0.82]
pvals /= pol_eff[0]
phi = 0.5*np.arctan2(U,Q)
dx = pvals*np.cos(phi)
dy = pvals*np.sin(phi)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
xs, ys = np.meshgrid(X,Y)
dx = dx[30:-30,230:-230]
dy = dy[30:-30,230:-230]
xs = xs[30:-30,230:-230]
ys = ys[30:-30,230:-230]
vectors = np.array([dx,dy])

def overplot(filein):
    img = mpimg.imread(filein)
    newim = np.zeros((470, 930, 3))
    newim = img[0:470,0:930,:]
    blast = fits.open(smooth250)
    wcs = WCS(blast[1].header)
    ctr_dec = Angle('-59d35m0s')
    ctr_ra = Angle('10h43m22s')
    fig = plt.figure(figsize = (14,13))
    f = aplpy.FITSFigure(smooth250, figure = fig)
    #f.set_yaxis_coord_type('latitude')
    f.tick_labels.set_xformat('dd.dd')
    #f.set_xaxis_coord_type('longitude')
    f.tick_labels.set_yformat('dd.dd')
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size=14)
    f.set_theme('publication')
    f.add_scalebar(15./60.) # arcmin
    f.scalebar.set_label('10 pc')
    f.scalebar.set_color('white')
    f.scalebar.set_corner('bottom right')
    f.scalebar.set_linewidth(2)
    f.scalebar.set_font_size(size = 'large')
    f.scalebar.set_label('10 pc')
    f.recenter(ctr_ra.deg, ctr_dec.deg, width = 1.5, height = 1.1)
    ax = plt.gca()
    ax.imshow(np.flipud(newim), origin = 'lower')
    vec = False
    streamlines = True
    if vec:
        v = '_vec'
        plot_vectors(ax, vectors, ys, xs, nskip = 10, alph = 0.3, col = 'red')
    else:
        v = ''
    if streamlines:
        sl = '_sl'
        plot_streams(ax, vectors, xs, ys, nskip = 25, alph = 0.4, col = 'yellow', vec = False)
    else:
        sl = ''
    plt.tight_layout()
    fileout = filein[filein.rfind('/') + 1:-4] + v + sl + '.png'
    plt.savefig(os.path.join(save_files_here, fileout), format='png', bbox_inches = 'tight')
    return

for f in files:
    overplot(f)

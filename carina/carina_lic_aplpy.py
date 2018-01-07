from getIQU import IQU
from subprocess import call
import sys
import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.visualization import AsinhStretch
from astropy.io import fits
import scipy.ndimage as ndimage
from skimage import filters
import aplpy
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, PowerStretch)
plt.ion()

if len(sys.argv) < 2:
    print "You must supply a band, 'e.g., carina_lic.py 250'"
    sys.exit()
else:
    band = sys.argv[1]

# Call this function with desired band as first argument, e.g.:
# python carina_lic.py 250
# lists over which to iterate:
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']
pol_eff = [0.81, 0.79, 0.82]

# define file paths
blastpol_dir = './carinaData/smooth/3.0_arcmin'
filename = glob.glob(blastpol_dir + '/carinaneb_' + band + '_smoothed_3.0_rl.fits')[0]

# load in I, Q, U for desired band
Ivals, Qvals, Uvals, wcs = IQU(band, filename)
wcs.wcs.crpix[0] -= 260
wcs.wcs.crpix[1] -= 30
I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

# Correct pvals as in Jamil's thesis, 5.7
#pvals[pvals > 0.5] = np.nan
pvals /= pol_eff[bands.index(band)]
phi = 0.5*np.arctan2(U,Q) + np.pi/4.
dx = np.cos(phi)
dy = np.sin(phi)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
xs, ys = np.meshgrid(X,Y)

"""
plt.figure()
nskip = 2
skip = (slice(None, None, nskip), slice(None, None, nskip))
#f = aplpy.FITSFigure(I, figsize = (10.24,7.68), dpi = 100)
ax = plt.gca()
ax.imshow(I, cmap = "gist_heat")
#f.tick_labels.set_font(size='small')
#f.show_colorscale(cmap='gist_heat')
# Add polarization vectors
ax.quiver(xs[skip],ys[skip],(dx/mag)[skip],(dy/mag)[skip], color = "white", angles = 'xy', units = 'xy', scale_units = 'xy', scale = 0.3)
#f.show_vectors(pvals, phi, color = 'white', rotate = 90., scale = 50, step = 10) 
ax.set_facecolor('black')
plt.tight_layout()
"""
xsize, ysize = len(X), len(Y)

vectors = np.array([dx,dy])
#white = np.random.rand(xsize, ysize)
#white = np.random.uniform(low = 0., high = 1., size = (xsize, ysize))
white = np.random.normal(0., 1., size = (xsize,ysize))
#white = scipy.ndimage.gaussian_filter(white, sigma)

with file('texture.dat', 'w') as outfile:
    for row in white:
        np.savetxt(outfile, row, newline = " ")
        outfile.write('\n')
with file('dx.dat', 'w') as outfile:
    for row in dx:
        np.savetxt(outfile, row, newline = " ")
        outfile.write('\n')
with file('dy.dat', 'w') as outfile:
    for row in dy:
        np.savetxt(outfile, row, newline = " ")
        outfile.write('\n')

# command = ["./carina_lic", str(xsize), str(ysize)]
# call(command)

lic = np.loadtxt("./lic.dat")
#np.save('lic.npy', lic)
lic = np.transpose(lic)
#lic = np.rot90(lic)
#lic = np.load("lic.npy")
#lic += np.abs(np.nanmin(lic))
#lic[lic > 3*np.nanstd(lic)] *= 100*lic[lic > 3*np.nanstd(lic)]
mult = lic * I

#stretch = SqrtStretch()
#lic_stretch = stretch(lic)
hdu1 = fits.PrimaryHDU(data=lic, header=wcs.to_header())
f1 = aplpy.FITSFigure(hdu1)
f1.set_theme('publication')
f1.show_colorscale(stretch='linear', cmap='inferno', smooth = None, kernel='gauss',
                                             aspect='equal', interpolation='hamming')
ax = plt.gca()
ax.set_facecolor("k")
f1.axis_labels.set_font(size=16)
f1.tick_labels.set_font(size=14)
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
plt.tight_layout()

hdu2 = fits.PrimaryHDU(data=I, header=wcs.to_header())
f2 = aplpy.FITSFigure(hdu2)
f2.set_theme('publication')
#ax = plt.gca()
ax.set_facecolor("k")
f2.show_colorscale(cmap = 'inferno')
f2.axis_labels.set_font(size=16)
f2.tick_labels.set_font(size = 14)
#  scalebar
f2.add_scalebar(30/60.) # arcmin
f2.scalebar.set_label('0.5 deg')
f2.scalebar.set_color('white')
f2.scalebar.set_corner('bottom right')
f2.scalebar.set_label('0.5 deg')
f2.scalebar.set_linewidth(2)
f2.scalebar.set_font_size(size = 'large')
f2.add_grid()
f2.grid.set_color('yellow')
f2.grid.set_alpha(0.2)
norm = ImageNormalize(lic, interval=MinMaxInterval(), stretch=PowerStretch(1.2))
plt.imshow(lic, alpha = 0.3, origin = 'lower', norm = norm, interpolation = 'hamming', cmap = 'gray')
plt.savefig('./lic_overplot.png', dpi = 100, bbox_inches = 'tight')
plt.tight_layout()

hdu3 = fits.PrimaryHDU(data=mult, header=wcs.to_header())
f3 = aplpy.FITSFigure(hdu3)
#plt.imshow(mult, origin = 'lower', interpolation = 'gaussian', cmap = 'inferno')
f3.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f3.show_colorscale(vmin=0., vmax=0.005, stretch='arcsinh', cmap='inferno',\
             smooth = None, kernel='gauss', aspect='equal', interpolation='hamming')
#f3.show_colorscale(vmin=0., vmax=0.00015, stretch='arcsinh', cmap='inferno',\
#             smooth = None, kernel='gauss', aspect='equal', interpolation='hamming')
f3.axis_labels.set_font(size=16)
f3.tick_labels.set_font(size = 14)
#  scalebar
f3.add_scalebar(30/60.) # arcmin
f3.scalebar.set_label('0.5 deg')
f3.scalebar.set_color('white')
f3.scalebar.set_corner('bottom right')
f3.scalebar.set_label('0.5 deg')
f3.scalebar.set_linewidth(2)
f3.scalebar.set_font_size(size = 'large')
f3.add_grid()
f3.grid.set_color('yellow')
f3.grid.set_alpha(0.2)
plt.savefig('./mult.png', dpi = 100, bbox_inches = 'tight')
plt.tight_layout()

#im = convolve_fft(lic, I, fft_pad=True, psf_pad=True)
#im = ndimage.sobel(mult)
#plt.figure()
#plt.imshow(im, origin = 'lower', cmap = 'inferno')


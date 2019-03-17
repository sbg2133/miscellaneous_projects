import sys, os
from getIQU import IQU
from subprocess import call
import sys, os
import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.visualization import AsinhStretch
from astropy.io import fits
from astropy.wcs import WCS
import scipy.ndimage as ndimage
from skimage import filters
import aplpy
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, PowerStretch)
plt.ion()
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"

stokes = ['I', 'Q', 'U']
planck_dir = './carinaData/planckData'
filename = os.path.join(planck_dir, "planck_353_carinaneb_pol.fits")
hdulist = fits.open(filename)
for s, param in enumerate(stokes):
    if (stokes[s] == 'I'):
        Ivals = hdulist[s+1].data
        Ivals[Ivals == 0.0] = np.nan
        #wcs = WCS(hdulist[s+1].header)
for s, param in enumerate(stokes):
    if (stokes[s] == 'Q'):
        Qvals = hdulist[s+1].data
        Qvals[Qvals == 0.0] = np.nan
        #wcs = WCS(hdulist[s+1].header)
for s, param in enumerate(stokes):
    if (stokes[s] == 'U'):
        Uvals = hdulist[s+1].data
        Uvals[Uvals == 0.0] = np.nan
        wcs = WCS(hdulist[s+1].header)

wcs.wcs.crpix[0] -= 260
wcs.wcs.crpix[1] -= 30
I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

phi = 0.5*np.arctan2(U,Q)
dx = pvals*np.cos(phi)
dy = pvals*np.sin(phi)
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
sigma = 1.2
white = ndimage.gaussian_filter(white, sigma)

with file('texture_planck.dat', 'w') as outfile:
    for row in white:
        np.savetxt(outfile, row, newline = " ")
        outfile.write('\n')
with file('dx_planck.dat', 'w') as outfile:
    for row in dx:
        np.savetxt(outfile, row, newline = " ")
        outfile.write('\n')
with file('dy_planck.dat', 'w') as outfile:
    for row in dy:
        np.savetxt(outfile, row, newline = " ")
        outfile.write('\n')

#command = ["./planck_lic", str(xsize), str(ysize)]
#call(command)

lic = np.loadtxt("./lic_planck.dat")
lic = np.transpose(lic)
#lic += np.abs(np.nanmin(lic))
mult = lic * I

hdu2 = fits.PrimaryHDU(data=np.zeros_like(I), header=wcs.to_header())
f2 = aplpy.FITSFigure(hdu2, figsize = (10,10))
f2.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f2.add_scalebar(15/60.) # arcmin
f2.scalebar.set_color('white')
f2.scalebar.set_corner('bottom right')
f2.scalebar.set_label('10 pc')
f2.tick_labels.set_yformat('dd.dd')
f2.tick_labels.set_xformat('dd.dd')
f2.axis_labels.set_font(size=16)
f2.tick_labels.set_font(size=14)
plt.imshow(I, origin = 'lower', cmap = "inferno", alpha = 1)
#plt.imshow(lic, vmin = -0.07, vmax = 0.3, origin = 'lower', cmap = "gray", alpha = 0.4, interpolation = "bilinear")
plt.imshow(lic, vmin = -0.05, vmax = 0.2, origin = 'lower', cmap = "gray", alpha = 0.5, interpolation = "bilinear")
plt.tight_layout()
#plt.savefig(os.path.join(save_files_here, 'planck_han_51.eps'), format='eps', dpi=1000, bbox_inches = 'tight')
plt.savefig(os.path.join(save_files_here, 'planck_han_51.png'), format='png', bbox_inches = 'tight')

hdu3 = fits.PrimaryHDU(data=np.zeros_like(I), header=wcs.to_header())
f3 = aplpy.FITSFigure(hdu3, figsize = (10,10))
f3.set_theme('publication')
#  scalebar
ax = plt.gca()
ax.set_facecolor("k")
f3.add_scalebar(15/60.) # arcmin
f3.scalebar.set_color('white')
f3.scalebar.set_corner('bottom right')
f3.scalebar.set_label('10 pc')
f3.tick_labels.set_yformat('dd.dd')
f3.tick_labels.set_xformat('dd.dd')
f3.axis_labels.set_font(size=16)
f3.tick_labels.set_font(size=14)
vmin = [-0.00007, 0.4, 0.5]
vmax = [0.02, 0.4, 0.5]
plt.imshow(mult, origin = 'lower',\
      cmap = "inferno", vmin = vmin[0],\
     vmax = vmax[0], interpolation = 'bilinear')
#plt.tight_layout()
#plt.savefig(os.path.join(save_files_here, 'planck2_han51.eps'), format='eps', dpi=500, bbox_inches = 'tight')
plt.savefig(os.path.join(save_files_here, 'planck2_han51.png'), bbox_inches = 'tight')
"""
#stretch = SqrtStretch()
#lic_stretch = stretch(lic)
hdu1 = fits.PrimaryHDU(data=lic, header=wcs.to_header())
f1 = aplpy.FITSFigure(hdu1)
f1.set_theme('publication')
f1.show_colorscale(stretch='linear', cmap='inferno', smooth = None, kernel='gauss', aspect='equal', interpolation='hamming')
ax = plt.gca()
ax.set_facecolor("k")
f1.axis_labels.set_font(size=16)
f1.tick_labels.set_font(size=14)
#  scalebar
f1.add_scalebar(15/60.) # arcmin
f1.scalebar.set_label('10 pc')
f1.scalebar.set_color('white')
f1.scalebar.set_corner('bottom right')
f1.scalebar.set_label('10 pc')
f1.scalebar.set_linewidth(2)
f1.scalebar.set_font_size(size = 'large')
#f1.add_grid()
#f1.grid.set_color('yellow')
#f1.grid.set_alpha(0.3)
plt.tight_layout()

hdu2 = fits.PrimaryHDU(data=I, header=wcs.to_header())
f2 = aplpy.FITSFigure(hdu2)
f2.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f2.show_colorscale(cmap = 'inferno')
#f2.show_colorscale(stretch='linear', cmap='inferno', aspect='equal', interpolation='hamming', vmin=0., vmax=0.0018)
f2.axis_labels.set_font(size=16)
f2.tick_labels.set_font(size = 14)
#  scalebar
f2.add_scalebar(15/60.) # arcmin
f2.scalebar.set_label('10 pc')
f2.scalebar.set_color('white')
f2.scalebar.set_corner('bottom right')
f2.scalebar.set_label('10 pc')
f2.scalebar.set_linewidth(2)
f2.scalebar.set_font_size(size = 'large')
#f2.add_grid()
#f2.grid.set_color('yellow')
#f2.grid.set_alpha(0.2)
norm = ImageNormalize(lic, interval=MinMaxInterval(), stretch=PowerStretch(1.05))
plt.imshow(lic, alpha = 0.4, origin = 'lower', interpolation = 'hamming', cmap = 'gray')
#plt.savefig('./lic_overplot.png', dpi = 100, bbox_inches = 'tight')
plt.tight_layout()

hdu3 = fits.PrimaryHDU(data=mult, header=wcs.to_header())
f3 = aplpy.FITSFigure(hdu3)
#plt.imshow(mult, origin = 'lower', interpolation = 'gaussian', cmap = 'inferno')
f3.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
#f3.show_colorscale(cmap = 'inferno')
f3.show_colorscale(vmin = 0., vmax = 0.6, stretch='linear', cmap='inferno',smooth = None, kernel='gauss', aspect='equal', interpolation='hamming')
#f3.show_colorscale(vmin=0., vmax=100, cmap='inferno',smooth = None, kernel='gauss', aspect='equal', interpolation='hamming')
f3.axis_labels.set_font(size=16)
f3.tick_labels.set_font(size = 14)
#  scalebar
f3.add_scalebar(15/60.) # arcmin
f3.scalebar.set_label('10 pc')
f3.scalebar.set_color('white')
f3.scalebar.set_corner('bottom right')
f3.scalebar.set_label('10 pc')
f3.scalebar.set_linewidth(2)
f3.scalebar.set_font_size(size = 'large')
#f3.add_grid()
#f3.grid.set_color('yellow')
#f3.grid.set_alpha(0.2)
#fits.writeto('./mult2.fits', mult, header=wcs.to_header())
#plt.savefig('./mult.png', dpi = 100, bbox_inches = 'tight')
plt.tight_layout()

#im = convolve_fft(lic, I, fft_pad=True, psf_pad=True)
#im = ndimage.sobel(mult)
#plt.figure()
#plt.imshow(im, origin = 'lower', cmap = 'inferno')
"""

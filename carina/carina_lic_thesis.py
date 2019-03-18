from getIQU import IQU
from subprocess import call
import sys,os
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
import scipy.ndimage
from makePretty import pretty
import aplpy
from skimage import filters
#matplotlib.rcParams.update({'font.size':16})
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"
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
blastpol_dir = './carinaData/smooth/5.0_arcmin_no_kernnorm'
filename = glob.glob(blastpol_dir + '/carinaneb_' + band + '_smoothed_5.0_rl_cal.fits')[0]
#blastpol_dir = './carinaData/smooth/3.0_arcmin'
#filename = glob.glob(blastpol_dir + '/carinaneb_' + band + '_smoothed_3.0_rl.fits')[0]

# load in I, Q, U for desired band
#Ivals, Qvals, Uvals, __, wcs = IQU(filename)
Ivals, Qvals, Uvals, pol_data, wcs = IQU(filename, do_cov = True)
phi2 = np.deg2rad(pol_data[1][30:-30,260:-260])
p = pol_data[2][30:-30,260:-260]
p[p > 0.5] = np.nan
p[p < 0] = np.nan

I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

# Correct pvals as in Jamil's thesis, 5.7
pvals[pvals > 0.5] = np.nan
pvals[pvals < 0] = np.nan
pvals /= pol_eff[bands.index(band)]
phi = 0.5*np.arctan2(U,Q)
#dx = np.cos(phi)
#dy = np.sin(phi)
#dx = pvals*np.cos(phi)
#dy = pvals*np.sin(phi)
dx = np.cos(phi2)
dy = np.sin(phi2)
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
white = scipy.ndimage.gaussian_filter(white, sigma)

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

if len(sys.argv) > 2:
    print
    print "Doing LIC"
    print
    command = ["./carina_lic", str(xsize), str(ysize)]
    call(command)

lic = np.loadtxt("./lic.dat")
lic = np.transpose(lic)
lic -= np.nanmin(lic)
#np.save('lic.npy', lic)
#lic2 = np.load("lic.npy")
#lic = scipy.ndimage.gaussian_filter(lic, 1.01)
mult = lic * I
mult -= np.nanmin(lic)

"""
blur_size = 8
unsharp_strength = 0.8
blurred = filter.gaussian_filter(lic, blur_size)
highpass = lic - unsharp_strength*blurred
sharp = lic + highpass
lowpass = scipy.ndimage.gaussian_filter(lic, 5)
highpass = lic - lowpass
highpass += lic
"""
hdu2 = fits.PrimaryHDU(data=np.zeros_like(I), header=wcs.to_header())
f2 = aplpy.FITSFigure(hdu2, figsize = (10,10))
f2.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
f2.add_scalebar(15/60.) # arcmin
f2.scalebar.set_label('0.5 deg')
f2.scalebar.set_color('white')
f2.scalebar.set_corner('bottom right')
f2.scalebar.set_label('10 pc')
f2.tick_labels.set_yformat('dd.dd')
f2.tick_labels.set_xformat('dd.dd')
f2.axis_labels.set_font(size=16)
f2.tick_labels.set_font(size=14)
plt.imshow(I, origin = 'lower', cmap = "inferno", alpha = 1)
#plt.imshow(lic, vmin = -0.07, vmax = 0.3, origin = 'lower', cmap = "gray", alpha = 0.4, interpolation = "bilinear")
plt.imshow(lic, vmin = 0.25, vmax = 1.0, origin = 'lower', cmap = "gray", alpha = 0.5, interpolation = "bilinear")
plt.tight_layout()
#f2.savefig(os.path.join(save_files_here, 'lic_han_51.eps'), format='eps', dpi=1000, transparent = True)
#plt.savefig(os.path.join(save_files_here, 'lic_han_51.png'), format='png', bbox_inches = 'tight')

##################################################
# For 250 um: v = 1000
# For 350 um: v = 500
# For 500 um: v = 200

hdu3 = fits.PrimaryHDU(data=np.zeros_like(I), header=wcs.to_header())
f3 = aplpy.FITSFigure(hdu3, figsize = (10,10))
f3.set_theme('publication')
#  scalebar
ax = plt.gca()
ax.set_facecolor("k")
f3.add_scalebar(15/60.) # arcmin
f3.scalebar.set_label('0.5 deg')
f3.scalebar.set_color('white')
f3.scalebar.set_corner('bottom right')
f3.scalebar.set_label('10 pc')
f3.tick_labels.set_yformat('dd.dd')
f3.tick_labels.set_xformat('dd.dd')
f3.axis_labels.set_font(size=16)
f3.tick_labels.set_font(size=14)
vmin = [5, 5, 5]
vmax = [200, 200, 200]
plt.imshow(mult, origin = 'lower',\
      cmap = "inferno", vmin = vmin[bands.index(band)],\
     vmax = vmax[bands.index(band)], interpolation = 'bilinear')
#plt.tight_layout()
#plt.savefig(os.path.join(save_files_here, 'lic2_han51.png'), format='png', bbox_inches = 'tight')

"""
hdu4 = fits.PrimaryHDU(data=np.zeros_like(I), header=wcs.to_header())
f4 = aplpy.FITSFigure(hdu4, figsize = (10,10))
f4.set_theme('publication')
#  scalebar
ax = plt.gca()
ax.set_facecolor("k")
f4.add_scalebar(15/60.) # arcmin
f4.scalebar.set_color('white')
f4.scalebar.set_corner('bottom right')
f4.scalebar.set_label('10 pc')
f4.tick_labels.set_yformat('dd.dd')
f4.tick_labels.set_xformat('dd.dd')
f4.axis_labels.set_font(size=16)
f4.tick_labels.set_font(size=16)
vmin = [-0.000007, 0.4, 0.5]
vmax = [0.0001, 0.4, 0.5]
plt.imshow(mult, origin = 'lower',\
      cmap = "viridis", vmin = vmin[bands.index(band)],\
     vmax = vmax[bands.index(band)], interpolation = 'bilinear')
#plt.tight_layout()
#plt.savefig(os.path.join(save_files_here, 'lic2_viridis.png'), format='png', bbox_inches = 'tight')

"""

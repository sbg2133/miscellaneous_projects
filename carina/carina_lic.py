from getIQU import IQU
from subprocess import call
import sys
import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
import scipy.ndimage
from skimage import filters
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
Ivals, Qvals, Uvals, __, wcs = IQU(filename)

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
sigma = 1.02
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

command = ["./carina_lic", str(xsize), str(ysize)]
call(command)

lic = np.loadtxt("./lic.dat")
np.save('lic.npy', lic)
lic = np.transpose(lic)
lic2 = np.load("lic.npy")
lic2 += np.abs(np.nanmin(lic2))
lic2[lic2 > 3*np.nanstd(lic2)] *= 100*lic2[lic2 > 3*np.nanstd(lic2)]
mult =  lic2 * I

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
"""
fig1 = plt.figure(figsize=(10.24, 7.68), dpi = 100)
ax = fig1.add_subplot(1, 1, 1, projection=wcs)
ax.set_facecolor("k")
#plt.pcolor(X, Y, lic, cmap='inferno')
ax.imshow(lic, cmap = "inferno", interpolation = "hanning")
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel('RA', fontsize = 16, fontweight = 'bold')
ax.set_ylabel('DEC', fontsize = 16, fontweight = 'bold')
plt.tight_layout()
"""

fig2 = plt.figure()
ax = fig2.add_subplot(1, 1, 1, projection=wcs)
plt.imshow(I, cmap = "inferno", alpha = 0.9)
plt.imshow(lic, cmap = "gray", alpha = 0.4, interpolation = "bilinear")
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel('RA', fontsize = 18)
ax.set_ylabel('DEC', fontsize = 18)
#ax.set_facecolor("k")
#plt.tight_layout()

##################################################
# For 250 um: v = 1000
# For 350 um: v = 500
# For 500 um: v = 200

v = [0.25, 0.4, 0.5]
fig3 = plt.figure()
ax = fig3.add_subplot(1, 1, 1, projection=wcs)
plt.imshow(mult, cmap = "inferno", vmin = 0, vmax = v[bands.index(band)])
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel('RA', fontsize = 18)
ax.set_ylabel('DEC', fontsize = 18)
#ax.set_facecolor("k")
#plt.tight_layout()

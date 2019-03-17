from getIQU import IQU
from subprocess import call
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.wcs import WCS
import scipy.ndimage
from skimage import filters
plt.ion()

stokes = ['I', 'Q', 'U']

planck_dir = './carinaData/planckData'
filename = os.path.join(planck_dir, "planck_353_carinaneb_pol.fits")
hdulist = fits.open(filename)
for s, param in enumerate(stokes):
    if (stokes[s] == 'I'):
        I = hdulist[s+1].data
        I[I == 0.0] = np.nan
        #wcs = WCS(hdulist[s+1].header)
for s, param in enumerate(stokes):
    if (stokes[s] == 'Q'):
        Q = hdulist[s+1].data
        Q[Q == 0.0] = np.nan
        #wcs = WCS(hdulist[s+1].header)
for s, param in enumerate(stokes):
    if (stokes[s] == 'U'):
        U = hdulist[s+1].data
        U[U == 0.0] = np.nan
        wcs = WCS(hdulist[s+1].header)

I = I[30:-30,260:-260]
Q = Q[30:-30,260:-260]
U = U[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

# Correct pvals as in Jamil's thesis, 5.7
#pvals[pvals > 0.5] = np.nan
phi = 0.5*np.arctan2(U,Q)
dx = np.cos(phi)
dy = np.sin(phi)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
print "Y =", I.shape[1]
print "X =", I.shape[0]
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

command = ["./carina_lic", str(ysize), str(xsize)]
call(command)

lic = np.loadtxt("./lic.dat")
#np.save('lic.npy', lic)
lic = np.transpose(lic)
#lic = np.load("lic.npy")
lic += np.abs(np.nanmin(lic))
lic[lic > 3*np.nanstd(lic)] *= 100*lic[lic > 3*np.nanstd(lic)]
mult =  lic * I

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
fig1 = plt.figure(figsize=(10.24, 7.68), dpi = 100)
ax = fig1.add_subplot(1, 1, 1, projection=wcs)
ax.set_facecolor("k")
ax.imshow(lic, cmap = "inferno", interpolation = "hanning")
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel('RA', fontsize = 16, fontweight = 'bold')
ax.set_ylabel('DEC', fontsize = 16, fontweight = 'bold')
plt.tight_layout()

fig2 = plt.figure(figsize=(10.24, 7.68), dpi = 100)
ax = fig2.add_subplot(1, 1, 1, projection=wcs)
plt.imshow(I, cmap = "inferno", alpha = 1)
plt.imshow(lic, cmap = "gray", alpha = 0.30, interpolation = "hanning")
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel('RA', fontsize = 16, fontweight = 'bold')
ax.set_ylabel('DEC', fontsize = 16, fontweight = 'bold')
ax.set_facecolor("k")
plt.tight_layout()

fig3 = plt.figure(figsize=(10.24, 7.68), dpi = 100)
ax = fig3.add_subplot(1, 1, 1, projection=wcs)
plt.imshow(mult, cmap = "inferno", vmin = 0, vmax = 100.)
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel('RA', fontsize = 16, fontweight = 'bold')
ax.set_ylabel('DEC', fontsize = 16, fontweight = 'bold')
ax.set_facecolor("k")
plt.tight_layout()

from subprocess import call
import sys
import numpy as np
import glob
import matplotlib.pyplot as plt
import aplpy
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
import scipy.ndimage
from skimage import filter

# lists over which to iterate:
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']
pol_eff = [0.81, 0.79, 0.82]

# define file paths
blastpol_dir = '../'
filenames = glob.glob(blastpol_dir+'*p10_good_C_gls_map_cal.fits')

fwhm_orig = 36. # 250um fwhm, "
omega_orig = 1.13*fwhm_orig**2
smooth_fwhm = 4.5 * 60. # 4.5' in "
fwhm_norm = smooth_fwhm / fwhm_orig
sigma = fwhm_norm/2.3548
gauss_kernel = Gaussian2DKernel(sigma)
percent= 0.075 # empirically-tweaked

for b, filename in enumerate(filenames):
    # Read in the HDUlist from the FITS file:
    print(filename)
    hdulist = fits.open(filename)
    # Print some info on what's in it:
    hdulist.info()
    for s, param in enumerate(stokes):
        if (stokes[s] == 'I'):
            Ivals_raw = hdulist[s+1].data
            minval = np.nanmin(Ivals_raw)
            maxval = np.nanmax(Ivals_raw)
            # Subtracting mean flux from dark region
            far_region = Ivals_raw[120:350, 750:800]
            Ivals = Ivals_raw - np.mean(far_region)
            #Ivals[Ivals < np.std(far_region)] = np.nan
            Ivals = convolve(Ivals, gauss_kernel)
            Ivals[Ivals == 0.0] = np.nan
        if (stokes[s] == 'Q'):
            Qvals_raw = hdulist[s+1].data
            far_region = Qvals_raw[120:350, 750:800]
            Qvals = Qvals_raw - np.mean(far_region)
            #Qvals[Qvals < np.std(far_region)] = np.nan
            Qvals = convolve(Qvals, gauss_kernel)
            Qvals[Qvals == 0.0] = np.nan
        if (stokes[s] == 'U'):
            Uvals_raw = hdulist[s+1].data
            far_region = Uvals_raw[120:350, 750:800]
            Uvals = Uvals_raw - np.mean(far_region)
            #Uvals[Uvals < np.std(far_region)] = np.nan
            Uvals = convolve(Uvals, gauss_kernel)
            Uvals[Uvals == 0.0] = np.nan
I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
"""
I = scipy.ndimage.zoom(I, 2, order=0)
Q = scipy.ndimage.zoom(Q, 2, order=0)
U = scipy.ndimage.zoom(U, 2, order=0)
"""
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I

# Correct pvals as in Jamil's thesis, 5.7
#pvals[pvals > 0.5] = np.nan
pvals /= pol_eff[0]
phi = 0.5*np.arctan2(U,Q)
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
white = np.random.rand(xsize, ysize)
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
lic = np.transpose(lic)
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
plt.figure(figsize=(10.24, 7.68), dpi = 100)
ax = plt.gca()
ax.set_facecolor("k")
ax.imshow(lic, cmap = "inferno", interpolation = "hanning")
plt.tight_layout()

plt.figure(figsize=(10.24, 7.68), dpi = 100)
plt.imshow(I, cmap = "inferno", alpha = 1)
plt.imshow(lic, cmap = "gray", vmin = 0.35, alpha = 0.30, interpolation = "hanning")
ax = plt.gca()
ax.set_facecolor("k")
plt.tight_layout()

plt.figure(figsize=(10.24, 7.68), dpi = 100)
plt.imshow(mult, cmap = "inferno", vmin=0, vmax=2000)
ax = plt.gca()
ax.set_facecolor("k")
plt.tight_layout()



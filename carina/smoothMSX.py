import os
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.convolution import Gaussian2DKernel, convolve_fft, convolve
from skimage import restoration
import matplotlib.pyplot as plt
import blast.util

fwhm = 3.0
mapdir = "/home/wizwit/miscellaneous_projects/carina/carinaData"
#mapfile = os.path.join(mapdir, "carinaneb/carinaneb_good_250_p10_good_C_gls_map_cal.fits")

mapfile = os.path.join(mapdir, "msx/msxmapA.fits")
kernel_size = None

# open FITS file and get map size
hdu = fits.open(mapfile)[0]
if kernel_size is not None:
  nx = ny = kernel_size
else:
  # number of x pixels
  nx  = hdu.header["NAXIS1"]
  print "Number of X pixels:", nx
  # numer of y pixels
  ny  = hdu.header["NAXIS2"]
  print "Number of Y pixels:", ny
# distance between pixels
pix = hdu.header["CDELT2"] * 3600. # arcseconds
print "Delta Pixel (arcsec):", pix

h = hdu.header
w = wcs.WCS(h)
x = np.mgrid[1:h["NAXIS2"]+1, 1:h["NAXIS1"]+1][1]
y = np.mgrid[1:h["NAXIS2"]+1, 1:h["NAXIS1"]+1][0]
ra, dec = w.wcs_pix2world(x, y, 1)

# construct gaussian
#gauss = blast.util.get_kernel(fwhm, pix_size=pix, map_size=(nx, ny))
# convolve map by the kernel
#print("\t...convolving map")
#blast.util.smooth_map(hdu, gauss)
#w = wcs.WCS(hdu.header)
#h = hdu.header


#outfile = os.path.join(mapdir, "msx_test_smooth_" + str(fwhm) + "_arcmin.fits")
#hdus.writeto(outfile)

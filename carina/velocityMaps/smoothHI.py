import os
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.convolution import Gaussian2DKernel, convolve_fft, convolve
from skimage import restoration
import matplotlib.pyplot as plt
import blast.util

fwhm = 5.0
mapdir = "/home/wizwit/miscellaneous_projects/carina/carinaData"

mapfile = os.path.join(mapdir, "HI/carina.HI.fits")
kernel_size = None

hdu = fits.open(mapfile)
h = hdu[0].header
w = wcs.WCS(h)
if kernel_size is not None:
  nx = ny = kernel_size
else:
  # number of x pixels
  nx  = h["NAXIS1"]
  print "Number of X pixels:", nx
  # numer of y pixels
  ny  = h["NAXIS2"]
  print "Number of Y pixels:", ny
  nz  = h["NAXIS3"]
  print "Number of slices:", nz

# distance between pixels
pix = h["CDELT2"] * 3600. # arcseconds
print "Delta Pixel (arcsec):", pix
count = 0

x = np.mgrid[1:h["NAXIS2"]+1, 1:h["NAXIS1"]+1][1]
y = np.mgrid[1:h["NAXIS2"]+1, 1:h["NAXIS1"]+1][0]

# construct gaussian
gauss = blast.util.get_kernel(fwhm, pix_size=pix, map_size=(nx, ny))

smoothed_data = np.zeros_like(hdu[0].data[0])
while count < nz:
    print count
    # open FITS file and get map size
    new_hdu = fits.PrimaryHDU(hdu[0].data[0][count], header=w.to_header())
    # convolve map by the kernel
    print("\t...convolving map")
    blast.util.smooth_map(new_hdu, gauss)
    smoothed_data[count] = new_hdu.data
    count += 1

out_hdu = fits.PrimaryHDU(smoothed_data, header=w.to_header())
outfile = os.path.join('./', "HI_smooth_" + str(fwhm) + "_arcmin.fits")
out_hdu.writeto(outfile)

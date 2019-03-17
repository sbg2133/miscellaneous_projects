import os
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.convolution import Gaussian2DKernel, convolve_fft, convolve
from skimage import restoration
import matplotlib.pyplot as plt
import blast.util
plt.ion()

fwhm = 3.0
mapdir = "/home/wizwit/miscellaneous_projects/carina/carinaData/herschelData/SPIRE-P_HiRes"
blast_mapdir = "/home/wizwit/miscellaneous_projects/carina/carinaData"
blast_mapfile = os.path.join(blast_mapdir, "carinaneb/carinaneb_good_250_p10_good_C_gls_map_cal.fits")

mapfile = os.path.join(mapdir, "1342211615_hires_250_level25.fits")
kernel_size = None

# open FITS file and get map size
hdu = fits.open(mapfile)[1]
h = hdu.header
if kernel_size is not None:
  nx = ny = kernel_size
else:
  # number of x pixels
  nx  = h["naxis1"]
  print "Number of X pixels:", nx
  # numer of y pixels
  ny  = h["naxis2"]
  print "Number of Y pixels:", ny
# distance between pixels
pix = h["CDELT2"] * 3600. # arcseconds
print "Delta Pixel (arcsec):", pix

w = wcs.WCS(h)

# construct gaussian
#gauss = blast.util.get_kernel(fwhm, pix_size=pix, map_size=(nx, ny))
# convolve map by the kernel
#print("\t...convolving map")
#blast.util.smooth_map(hdu, gauss)

blast_hdu = fits.open(blast_mapfile)[1]
blast_h = blast_hdu.header
blast_w = wcs.WCS(blast_h)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=blast_w)
plt.imshow(hdu.data, cmap = "inferno", vmin = 0, vmax = 2000)
plt.imshow(blast_hdu.data, cmap = "viridis", alpha = 0.5)

#outfile = os.path.join(mapdir, "msx_test_smooth_" + str(fwhm) + "_arcmin.fits")
#hdus.writeto(outfile)

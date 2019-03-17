import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy
from astropy.wcs import WCS
import sys, os
from getIQU import IQU
from astropy import coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.interpolate import griddata
from makePretty import pretty
from mpl_toolkits.axes_grid1 import make_axes_locatable
import aplpy
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ion()

pol_eff = [0.81, 0.79, 0.82]
cen_coord = [160.84429, -59.582361]
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"
root_dir = '/home/wizwit/miscellaneous_projects/carina/carinaData'
b250 = 'smooth/5.0_arcmin_no_kernnorm/carinaneb_500_smoothed_5.0_rl_cal.fits'
b350 = 'smooth/5.0_arcmin_no_kernnorm/carinaneb_350_smoothed_5.0_rl_cal.fits'
b500 = 'smooth/5.0_arcmin_no_kernnorm/carinaneb_500_smoothed_5.0_rl_cal.fits'
b250 = os.path.join(root_dir, b250)
b350 = os.path.join(root_dir, b350)
b500 = os.path.join(root_dir, b500)

def getPsi(path_to_file, band_idx):
    I, Q, U, __, wcs = IQU(path_to_file)
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    pvals[pvals > 0.5] = np.nan
    pvals[pvals < 0] = np.nan
    pvals /= pol_eff[band_idx]
    psi = 0.5*np.arctan2(U,Q)
    return I, Q, U, wcs, pvals, psi

I250, Q250, U250, wcs, pvals, psi = getPsi(b250, 0)
I350, Q350, U350, wcs, pvals, psi = getPsi(b350, 1)
I500, Q500, U500, wcs, pvals, psi = getPsi(b500, 2)

"""
xs = np.arange(I250.shape[1])
ys = np.arange(I250.shape[0])
X,Y = np.meshgrid(xs, ys)
mask_ra, mask_dec = wcs.all_pix2world(X, Y, 0)
mask_gal_coords = SkyCoord(mask_ra*u.degree, mask_dec*u.degree, frame='fk5').galactic

wcs_gal = WCS(naxis=2)
wcs_gal.wcs.crpix = [563.0, 665.0]
wcs_gal.wcs.cdelt = np.array([-0.00277777798786, 0.00277777798786])
wcs_gal.wcs.crval = [287.41231676, -1.14784038739]
wcs_gal.wcs.ctype = ["GLON-SIN", "GLAT-SIN"]
hdu_gal = fits.PrimaryHDU(data=I250, header=wcs_gal.to_header())
fig = plt.figure()
fgal = aplpy.FITSFigure(hdu_gal, figure = fig)
#fgal.ticks.hide()
#fgal.tick_labels.hide()
#fgal.axis_labels.hide()
fgal.add_grid()
fgal.grid.set_color('black')
fgal.grid.show()
"""
# Intensity plot, 250 um
hdu = fits.PrimaryHDU(data=I250, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
#ax.set_facecolor("k")
#f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno')
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('red')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('black')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.set_axis_label_font(size=18)
#f.colorbar.set_axis_label_text(r'MJy sr$^{-1}$')
f.colorbar.set_font(size = 18)
f.colorbar.show()
plt.savefig(os.path.join(save_files_here, 'carina_I250.eps'), format='eps', dpi = 500, bbox_inches = 'tight')

# Intensity plot, 350 um
hdu = fits.PrimaryHDU(data=I350, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
#f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno')
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('red')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('black')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.set_axis_label_font(size=18)
#f.colorbar.set_axis_label_text(r'MJy str$^{-1}$')
f.colorbar.set_font(size = 18)
f.colorbar.show()
#plt.tight_layout()
plt.savefig(os.path.join(save_files_here, 'carina_I350.eps'), format='eps', dpi = 500, bbox_inches = 'tight')

# Intensity plot, 500 um
hdu = fits.PrimaryHDU(data=I500, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
#f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno')
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('red')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('black')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
#f.colorbar.set_axis_label_text(r'MJy str$^{-1}$')
f.colorbar.set_font(size = 18)
f.colorbar.set_axis_label_font(size=18)
f.colorbar.show()
#plt.tight_layout()
plt.savefig(os.path.join(save_files_here, 'carina_I500.eps'), format='eps', dpi = 500, bbox_inches = 'tight')

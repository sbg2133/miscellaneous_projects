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
from astropy.coordinates import Angle
from astropy import units as u
from scipy.interpolate import griddata
from add_regions import car_regions as reg
from makePretty import pretty
from mpl_toolkits.axes_grid1 import make_axes_locatable
import aplpy
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ion()
save = 1

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

# Intensity plot, 250 um
hdu = fits.PrimaryHDU(data=I250, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
#f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno')
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('white')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('white')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.set_axis_label_font(size=18)
#f.colorbar.set_axis_label_text(r'MJy sr$^{-1}$')
f.colorbar.set_font(size = 18)
f.colorbar.show()

reg(ax, f, wcs, c = 'lime', ls = '--', add_dots = False)

if save:
    plt.savefig(os.path.join(save_files_here, 'carina_I250.eps'), format='eps', dpi = 500, bbox_inches = 'tight')
"""
# Intensity plot, 350 um
hdu = fits.PrimaryHDU(data=I350, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
#f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno')
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('white')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('white')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
f.colorbar.set_axis_label_font(size=18)
#f.colorbar.set_axis_label_text(r'MJy str$^{-1}$')
f.colorbar.set_font(size = 18)
f.colorbar.show()
#plt.tight_layout()
if save:
    plt.savefig(os.path.join(save_files_here, 'carina_I350.eps'), format='eps', dpi = 500, bbox_inches = 'tight')

# Intensity plot, 500 um
hdu = fits.PrimaryHDU(data=I500, header=wcs.to_header())
f = aplpy.FITSFigure(hdu)
f.set_theme('publication')
ax = plt.gca()
ax.set_facecolor("k")
#f.recenter(cen_coord[0], cen_coord[1], width = 1.4, height = 1.25)
f.show_colorscale(cmap = 'inferno')
f.add_scalebar(15/60.) # arcmin
f.scalebar.set_color('white')
f.scalebar.set_corner('bottom right')
f.scalebar.set_label('10 pc')
f.scalebar.set_linewidth('2')
f.scalebar.set_font_size('16')
f.tick_labels.set_yformat('dd.dd')
f.tick_labels.set_xformat('dd.dd')
f.axis_labels.set_font(size=16)
f.tick_labels.set_font(size=16)
f.ticks.set_color('white')
f.ticks.set_linewidth('2')
f.add_colorbar()
f.colorbar.set_location('right')
#f.colorbar.set_axis_label_text(r'MJy str$^{-1}$')
f.colorbar.set_font(size = 18)
f.colorbar.set_axis_label_font(size=18)
f.colorbar.show()
#plt.tight_layout()
if save:
    plt.savefig(os.path.join(save_files_here, 'carina_I500.eps'), format='eps', dpi = 500, bbox_inches = 'tight')
"""
"""
fig = plt.figure()
ax = fig.add_subplot(111)

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel('RA (J2000)')
ax.set_ylabel('Dec (J2000)')

ax1 = fig.add_subplot(311, projection = wcs)
ax1.imshow(I250, origin = 'lower', cmap = 'inferno')
ax1.set_facecolor("k")
#pretty()
ax1.tick_params(color = 'white', top='on', bottom='on', left='on', right='on', width = 2)

ax2 = fig.add_subplot(312, projection = wcs)
ax2.imshow(I350, origin = 'lower', cmap = 'inferno')
#pretty()
ax2.set_facecolor("k")
ax2.tick_params(color = 'white', top='on', bottom='on', left='on', right='on', width = 2)

ax3 = fig.add_subplot(313, projection = wcs)
ax3.imshow(I500, origin = 'lower', cmap = 'inferno')
ax3.set_facecolor("k")
#pretty()
ax3.tick_params(color = 'white', top='on', bottom='on', left='on', right='on', width = 2)

fig.subplots_adjust(hspace=0)
"""

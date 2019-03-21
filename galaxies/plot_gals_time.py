import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sens_calc as sens
import sys, os
from astropy.io import fits
import aplpy
from astropy.wcs import WCS
plt.ion()

telescope = sys.argv[1]
rootdir = '.'
bands = np.array([250, 350, 500]) # microns
t_map = np.array([10., 5., 3., 1]) # hours
t_dust = 15.0 # K
mapsize = 0.5**2 # deg^2
blast06_beam_area = sens.arcsec2toDeg2(sens.A_beam_06) # deg^2
blastTNG_beam_area = sens.arcsec2toDeg2(sens.A_beam_tng) # deg^2
blast_beam_per_pix = blast06_beam_area / blastTNG_beam_area

iras_100_beam = sens.arcsec2toDeg2(2.0 * 60.0**2) # deg^2

sig_p = 0.01 # 1 percent

def galContour(obj_name, filepath, band_idx, IRAS = False):
    print filepath
    print "Band =", bands[band_idx]
    # Read in the HDUlist from the FITS file:
    hdulist = fits.open(filepath)
    # Print some info on what's in it:
    # Read in the data from the Image HDU:
    data = hdulist[0].data
    mapsize = 0.5**2 # deg^2

    # Minimum observable intensity (MJy/str) achievable by BLAST-TNG
    # by observing a map of map_size (deg^2), with 1-sigma error bars in
    # pol fraction, sig_p, in t hours."""
    contour_levels = sens.I_min(mapsize, sig_p, t_map, band_idx)
    print contour_levels
    Npix = data[data >= np.min(contour_levels)].size
    print "N pix:", Npix
    N_blast_beams = Npix * blast_beam_per_pix
    if IRAS:
        N_blast_beams = Npix * iras_100_beam / blastTNG_beam_area
        f_ratio = sens.f_nu_T(sens.c/(bands[band_idx] * 1.0e-6), t_dust) /\
          sens.f_nu_T(sens.c/(100.0e-6), t_dust)
        print "f(100)/f(" + str(bands[band_idx]) + ") =", f_ratio
        data *= f_ratio

    minval = np.nanmin(data)
    maxval = np.nanmax(data)
    sig = np.std(data)
    mean = np.mean(data)
 
    print "I min, I max =", maxval, minval, "MJy" 
    print "Number of BLAST beams =", N_blast_beams[band_idx]
    for i in range(len(t_map)):
        print "t (hr) =", t_map[i], ": I =",contour_levels[i], "MJy/str"
    print
    return contour_levels

def contourPlot(obj_name, filepath, band_idx, contour_levels):
    fig = plt.figure(figsize = (10.24, 7.68))
    hdulist = fits.open(filepath)
    data = hdulist[0].data
    wcs = WCS(hdulist[0].header)
    cen_coord = wcs.wcs.crval
    hdu = fits.PrimaryHDU(data, header=wcs.to_header())
    f = aplpy.FITSFigure(hdu, figure = fig)
    f.recenter(cen_coord[0], cen_coord[1], width = 40.0, height = 40.0)
    f.show_contour(filepath, cmap = matplotlib.cm.viridis, filled = True, levels = contour_levels, alpha = 0.5)
    f.set_theme('publication')
    ax = plt.gca()
    ax.set_facecolor("k")
    f.show_colorscale(cmap = 'gray')
    f.frame.set_linewidth(1)  # points
    f.frame.set_color('black')
    f.set_yaxis_coord_type('latitude')
    f.tick_labels.set_xformat('dd.dd')
    f.set_xaxis_coord_type('longitude')
    f.tick_labels.set_yformat('dd.dd')
    f.axis_labels.set_font(size = 16)
    f.tick_labels.set_font(size = 14)
    f.add_colorbar()
    f.colorbar.set_location('bottom')
    f.colorbar.set_axis_label_font(size = 'large')
    f.colorbar.set_axis_label_text('MJy sr^-1')
    f.colorbar.set_axis_label_rotation(90)
    f.colorbar.show()
    f.set_title(obj_name + ", " + str(bands[band_idx]) + r" $\mu$m Intensity", size = 16)
    #plt.tight_layout()
    return f

if telescope == 'blast':
    for i in range(len(bands)):
        #levels = galContour('NGC 1808', rootdir + '/ngc1808_' + str(bands[i]) + '.fits', i)
        #contourPlot('NGC 1808', rootdir + '/ngc1808_' + str(bands[i]) + '.fits', i, levels)
        levels = galContour('NGC 1566', rootdir + '/ngc1566_' + str(bands[i]) + '.fits', i)
        f = contourPlot('NGC 1566', rootdir + '/ngc1808_' + str(bands[i]) + '.fits', i, levels)

if telescope == 'iras':
    for i in range(len(bands)):
        galContour('M 83', rootdir + '/m83_iras_100um.fits', i, IRAS = True)
        galContour('NGC 4945', rootdir + '/ngc4945_iras_100um.fits', i, IRAS = True)
        galContour('NGC 1808', rootdir + '/ngc1808_iras_100um.fits', i, IRAS = True)
        galContour('NGC 6744', rootdir + '/ngc6744_iras_100um.fits', i, IRAS = True)
        galContour('NGC 5068', rootdir + '/ngc5068_iras_100um.fits', i, IRAS = True)
        galContour('NGC 2835', rootdir + '/ngc2835_iras_100um.fits', i, IRAS = True)

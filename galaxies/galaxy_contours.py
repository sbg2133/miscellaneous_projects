import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import sens_calc as sens
import sys, os
from astropy.io import fits
from astropy.wcs import WCS
from makePretty import pretty
import const as const

h = const.MKS_H
c = const.MKS_C
k = const.MKS_K

plt.ion()

telescope = sys.argv[1]
rootdir = '.'
bands = np.array([250, 350, 500]) # microns
#t_map = np.array([10., 5., 3., 1]) # hours
t_map = np.arange(1.0, 15.0) # hours
if telescope == 'iras':
    t_map = np.arange(1.0, 6.0) # hours
t_map = t_map[::-1]
t_dust = 10.0 # K
mapsize = 0.5**2 # deg^2
blast06_beam_area = sens.arcsec2toDeg2(sens.A_beam_06) # deg^2
blastTNG_beam_area = sens.arcsec2toDeg2(sens.A_beam_tng) # deg^2
blast_beam_per_pix = blast06_beam_area / blastTNG_beam_area

iras_100_beam = sens.arcsec2toDeg2(2.0 * 60.0**2) # deg^2

sig_p = 0.005 # 0.5 percent

def gal_contour(obj_name, filepath, band_idx, IRAS = False):
    print filepath
    print "Band =", bands[band_idx]
    # Read in the HDUlist from the FITS file:
    hdulist = fits.open(filepath)
    # Print some info on what's in it:
    # Read in the data from the Image HDU:
    if IRAS:
        data = hdulist[0].data
        minval = np.nanmin(data)
        maxval = np.nanmax(data)
    else:
        data = hdulist[0].data
        minval = np.nanmin(data)
        maxval = np.nanmax(data)
        xmin = np.where(data == maxval)[0][0] - 15
        ymin = np.where(data == maxval)[1][0] - 15
        xmax = np.where(data == maxval)[0][0] + 15
        ymax = np.where(data == maxval)[1][0] + 15
        data = data[xmin:xmax,ymin:ymax]
    mapsize = 0.5**2 # deg^2
    # Minimum observable intensity (MJy/str) achievable by BLAST-TNG
    # by observing a map of map_size (deg^2), with 1-sigma error bars in
    # pol fraction, sig_p, in t hours."""
    contour_levels = sens.I_min(mapsize, sig_p, t_map, band_idx)
    if IRAS:
        f_ratio = sens.Bmod(sens.c/(bands[band_idx] * 1.0e-6), t_dust) /\
          sens.Bmod(sens.c/(100.0e-6), t_dust)
        print 'f' + '(' + str(bands[band_idx]) + ')' + '/f(100 um) =', f_ratio
        scaled_data = data*f_ratio
        Npix = scaled_data[scaled_data >= np.min(contour_levels)].size
    else:
        Npix = data[data >= np.min(contour_levels)].size
    print "N pix:", Npix
    N_blast_beams = Npix * blast_beam_per_pix
    print "I min, I max =", maxval, minval, "MJy/sr" 
    print "Number of BLAST beams =", N_blast_beams[band_idx]
    for i in range(len(t_map)):
        print "t (hr) =", t_map[i], ": I =",contour_levels[i], "MJy/sr"
    print
    return data, contour_levels

def plotContours(obj_name, band_idx, data, contour_levels, IRAS = True):
    minval = np.nanmin(data)
    maxval = np.nanmax(data)
    if IRAS:
        deg_per_pix = 0.025 # deg per pixel
        arcmin_per_pix = deg_per_pix * 60.
        x = np.arange(-data.shape[0]/2, data.shape[0]/2)*arcmin_per_pix
    else:
        deg_per_pix = -0.0025
        arcmin_per_pix = deg_per_pix * 60.
        x = np.arange(-data.shape[0]/2, data.shape[0]/2)*arcmin_per_pix
        #hdulist = fits.open(filepath)
        #wcs = WCS(hdulist[0].header)
        #xs = np.arange(data.shape[1])
        #ys = np.arange(data.shape[0])
        #X, Y = np.meshgrid(xs, ys)
        #ra, dec = wcs.all_pix2world(X, Y, 0)
        #radec = SkyCoord(ra*u.degree, dec*u.degree, frame='fk5')
    y = x
    X,Y = np.meshgrid(x,y)
    fig = plt.figure(figsize=(10.24, 7.68), dpi=100)
    toextend = 'neither'
    colormap = 'viridis'
    #plt.imshow(data, origin='lower', interpolation='nearest',\
                 #cmap='gray', vmin=minval, vmax=maxval)
    plt.pcolor(X, Y, data, cmap='gray')
    cb_t = plt.colorbar(orientation='horizontal', extend=toextend)
    cb_t.ax.tick_params(labelsize = 14)
    cb_t.set_label(r"MJy sr$^{-1}$", size = 16)
    cb_t.set_ticklabels(t_map.astype('int'))
    plt.ylabel('arcmin', fontsize = 16)
    ax = plt.gca()
    ax.set_facecolor('black')
    cax = plt.contourf(X, Y, data, levels=(contour_levels),\
                 alpha = 0.5, cmap=colormap)
    #cb_t = plt.colorbar(cax, orientation='vertical',\
                 #spacing = 'proportional', extend=toextend)
    #cb_t.ax.tick_params(labelsize = 14)
    #cb_t.set_label("Observing time (hr)", size = 16)
    #cb_t.set_ticklabels(t_map.astype('int'))
    cb_I = plt.colorbar(cax, ticks = contour_levels,\
             orientation='vertical', extend=toextend)
    cb_I.ax.tick_params(labelsize = 14)
    labels = []
    for i in range(len(contour_levels)):
        labels.append(str(contour_levels[i].astype('int')) + '/' +\
           str(t_map[i].astype('int')))
    cb_I.ax.set_yticklabels(labels)
    cb_I.set_label(r"MJy sr$^{-1}$/min obs. time (hr)", size = 16)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
    plt.title(obj_name + ", " + str(bands[band_idx]) + r" $\mu$m Intensity", size = 16)
    plt.tight_layout()
    plt.savefig("./" + obj_name + "_" + str(bands[band_idx]) +'um.png', bbox_inches = 'tight')
    return

if telescope == 'blast':
    for i in range(len(bands)):
        data, contour_levels = gal_contour('NGC 1808', rootdir + '/ngc1808_' + str(bands[i]) + '.fits', i)
        plotContours('NGC 1808', i, data, contour_levels, IRAS = False)
        data, contour_levels = gal_contour('NGC 1566', rootdir + '/ngc1566_' + str(bands[i]) + '.fits', i)
        plotContours('NGC 1566', i, data, contour_levels, IRAS = False)

if telescope == 'iras':
    for i in range(len(bands)):
        obj = 'M 83'
        data, contour_levels = gal_contour(obj, rootdir + '/m83_iras_100um.fits', i, IRAS = True)
        plotContours(obj, i, data, contour_levels)
        data, contour_levels = gal_contour('NGC 4945', rootdir + '/ngc4945_iras_100um.fits', i, IRAS = True)
        plotContours('NGC 4945', i, data, contour_levels)
        #data, contour_levels = gal_contour('NGC 1808', rootdir + '/ngc1808_iras_100um.fits', i, IRAS = True)
        data, contour_levels = gal_contour('NGC 6744', rootdir + '/ngc6744_iras_100um.fits', i, IRAS = True)
        plotContours('NGC 6744', i, data, contour_levels)
        data, contour_levels = gal_contour('NGC 5068', rootdir + '/ngc5068_iras_100um.fits', i, IRAS = True)
        plotContours('NGC 5068', i, data, contour_levels)
        data, contour_levels = gal_contour('NGC 2835', rootdir + '/ngc2835_iras_100um.fits', i, IRAS = True)
        plotContours('NGC 2835', i, data, contour_levels)

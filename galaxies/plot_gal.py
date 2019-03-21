import numpy as np
import matplotlib.pyplot as plt
import sens_calc as sens
import sys, os
from astropy.io import fits

rootdir = '.'
bands = np.array([250., 350., 500.])

def gal_contour(filepath, band_idx):
# file_path = absolute path to fits file

    f_ratio = sens.f_nu_T(sens.c/(100. * 1.0e-6), 20.)/\
        sens.f_nu_T(sens.c/(bands[band_idx] * 1.0e-6), 20.)
    print filepath
    print "Band =", bands[band_idx]
    # Read in the HDUlist from the FITS file:
    hdulist = fits.open(filepath)
    # Print some info on what's in it:
    # Read in the data from the Image HDU:
    datamap = hdulist[0].data
    datamap = datamap/f_ratio
    mapsize = 0.5**2 # deg^2
    minval = np.nanmin(datamap) # the minimum not including NaN pixels
    maxval = np.nanmax(datamap)
    iras_beam_deg = (0.25/len(datamap[0]))**2 *np.pi
    blast_beam_per_pix = iras_beam_deg/sens.A_beam_deg
    sig = np.std(datamap)
    mean = np.mean(datamap)
    contour_levels = np.array([mean + 1.5*sig, mean + 3*sig, mean + 5*sig])
    Npix = datamap[datamap >= mean].size
    N_blast_beams = Npix * blast_beam_per_pix
    t = sens.obs_time(mapsize, 0.01, contour_levels, band_idx)
    t = np.round(t, 2)
    print "Number of BLAST beams =", N_blast_beams[band_idx]
    for i in range(len(t)):
        print "t (hr) =", t[i], ": I =",contour_levels[i], "MJy/str"
    print
    fig = plt.figure(figsize=(10.24, 7.68), dpi=100)
    toextend = 'neither'
    colourmap = 'gist_heat'
    plt.imshow(datamap, origin='lower', interpolation='nearest',
        cmap='gist_gray', vmin=minval, vmax=maxval)
    
    # cb_im = plt.colorbar(orientation='horizontal', extend=toextend)
    # get current axis object:
    ax = plt.gca()
    
    # Since white is part of the colour map, have a different
    # background colour for NaN pixels (pixels where there is no data):
    ax.set_axis_bgcolor('gray')
    
    plt.contourf(datamap, levels=(contour_levels), cmap=colourmap, alpha=1)
    
    cb = plt.colorbar(orientation='horizontal', extend=toextend)
    
    # cb.set_label('BLASTPol Carina Nebula '+bands[b]+' Micron Stokes $'+
        # stokes[s]+'$ Intensity [MJy/sr]', size=14)
    
    # Save the figure:
    #plt.savefig('./blastpol_carinaneb_'+bands[b]+'um'+'_'+stokes[s]+'.png')
    
    #plt.show()
    return datamap

for i in range(len(bands)):
    gal_contour(rootdir + '/m83_iras_100um.fits', i)


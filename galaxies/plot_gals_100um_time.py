import numpy as np
import matplotlib.pyplot as plt
import sens_calc as sens
from astropy.io import fits

bands = np.array([250, 350, 500])
t = np.array([10., 5., 3., 1])

rootdir = "."

def gal_contour(obj_name, filepath, band_idx):
    f_ratio = sens.f_nu_T(sens.c/(100. * 1.0e-6), 15.) /\
              sens.f_nu_T(sens.c/(bands[band_idx] * 1.0e-6), 15.)
    print filepath
    print "Band =", bands[band_idx]
    print "f(100)/f(" + str(bands[band_idx]) + ") =", f_ratio
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

    contour_levels = sens.I_min(mapsize, 0.01, t, band_idx)
    Npix = datamap[datamap >= np.min(contour_levels)].size
    N_blast_beams = Npix * blast_beam_per_pix
    print "Max I, Min I =", maxval, minval
    print "Number of BLAST beams =", N_blast_beams[band_idx]
    for i in range(len(t)):
        print "t (hr) =", t[i], ": I =",contour_levels[i], "MJy/str"
    print
    fig = plt.figure(figsize=(10.24, 7.68), dpi=100)
    toextend = 'neither'
    colourmap = 'gist_heat'
    plt.imshow(datamap, origin='lower', interpolation='nearest',\
              cmap='gist_gray', vmin=minval, vmax=maxval)

    # cb_im = plt.colorbar(orientation='horizontal', extend=toextend)
    # get current axis object:
    ax = plt.gca()

    ax.set_axis_bgcolor('gray')

    plt.contourf(datamap, levels=(contour_levels), cmap=colourmap, alpha=1)

    cb_t = plt.colorbar(orientation='horizontal', extend=toextend, shrink = 0.7)
    cb_t.set_label("Observing time (hr)", size = 14)
    cb_t.set_ticklabels(t)
    cb_I = plt.colorbar(orientation='vertical', extend=toextend)
    cb_I.set_label("Intensity (Mjy/sr)", size = 14)
    cb_I.set_ticklabels(np.round(contour_levels,1))

    plt.title(obj_name + ", " + str(bands[band_idx]) + r" $\mu$m Intensity", size = 14)

    plt.savefig("./" + obj_name + "_" + str(bands[band_idx]) +'um.png', bbox_inches = 'tight')

    #plt.show()
    return datamap

for i in range(len(bands)):
    gal_contour('M 83', rootdir + '/m83_iras_100um.fits', i)
    gal_contour('NGC 4945', rootdir + '/ngc4945_iras_100um.fits', i)
    gal_contour('NGC 1808', rootdir + '/ngc1808_iras_100um.fits', i)
    gal_contour('NGC 6744', rootdir + '/ngc6744_iras_100um.fits', i)
    gal_contour('NGC 5068', rootdir + '/ngc5068_iras_100um.fits', i)
    gal_contour('NGC 2835', rootdir + '/ngc2835_iras_100um.fits', i)


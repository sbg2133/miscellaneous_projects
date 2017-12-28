import numpy as np
import matplotlib.pyplot as plt
from getIQU import IQU
from streamLines import plot_streams

##################################################################
# run like: plotIntensity(['250'], streamlines = True, nskip = 12)
##################################################################
plt.ion()
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']
pol_eff = [0.81, 0.79, 0.82]
all_bands = ['250', '350', '500']

I250, Q250, U250, wcs = IQU('250', './carinaData/carinaneb/carinaneb_good_250_p10_good_C_gls_map_cal.fits')
I350, Q350, U350, wcs = IQU('350', './carinaData/carinaneb/carinaneb_good_350_p10_good_C_gls_map_cal.fits')
I500, Q500, U500, wcs = IQU('500', './carinaData/carinaneb/carinaneb_good_500_p10_good_C_gls_map_cal.fits')

Is = [I250, I350, I500]
Qs = [Q250, Q350, Q500]
Us = [U250, U350, U500]
"""
I250[np.isnan(I350)] = 0.
I250[I250 < 0.] = 0.
I350[np.isnan(I350)] = 0.
I350[I350 < 0.] = 0.
I500[np.isnan(I500)] = 0.
I500[I500 < 0.] = 0.
I250[I250 < 5*np.std(I250)] = 0.
I350[I350 < 5*np.std(I350)] = 0.
I500[I500 < 5*np.std(I500)] = 0.
"""
def plotIntensity(bands_list, streamlines = False, nskip = 10):
    f, ax = plt.subplots(len(bands_list), 1, sharex = True, dpi = 100., squeeze=False, subplot_kw = {'projection': wcs})
    f.text(0.5, 0.04, 'RA', ha='center', fontsize = 12)
    f.text(0.04, 0.5, 'DEC', va='center', rotation='vertical', fontsize = 12)
    #f.tight_layout()
    [ax[bands.index(band), 0].imshow(Is[bands.index(band)], cmap = "inferno") for band in bands_list]
    #for i in range(len(ax)):
        #ax[i, 0].set_facecolor("k")
    if streamlines:
        [putStreamlines(ax[bands.index(band), 0], bands.index(band), nskip) for band in bands_list]
    plt.tight_layout()
    return

def putStreamlines(ax, band_idx, nskip):
    """
    I = Is[band_idx][30:-30,260:-260]
    Q = Qs[band_idx][30:-30,260:-260]
    U = Us[band_idx][30:-30,260:-260]
    """
    I = Is[band_idx]
    Q = Qs[band_idx]
    U = Us[band_idx]
    Pvals = np.sqrt(Q**2 + U**2)
    pvals = Pvals/I
    # Correct pvals as in Jamil's thesis, 5.7
    #pvals[pvals > 0.5] = np.nan
    pvals /= pol_eff[band_idx]
    phi = 0.5*np.arctan2(U,Q) 
    dx = np.cos(phi)
    dy = np.sin(phi)
    mag = np.sqrt(dx**2 + dy**2)
    X = np.linspace(0, I.shape[1], I.shape[1])
    Y = np.linspace(0, I.shape[0], I.shape[0])
    xs, ys = np.meshgrid(X,Y)
    xsize, ysize = len(X), len(Y)
    vectors = np.array([dx,dy])
    plot_streams(ax, vectors, xs, ys, nskip)
    return

#plot_streams(vectors, xs, ys, nskip = 10, vec = True)

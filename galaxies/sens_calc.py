import numpy as np
import const as const

A_prim_06 = np.pi * 0.9**2 # m^2
Ndet_06 = np.array([139, 88, 43])
fwhm_06 = np.array([36.,42.,60.])
A_beam_06 = 1.13 * fwhm_06**2 # arcsec^2

# TNG primary mirror diameter = 2.5 m
A_prim_tng = np.pi * (2.5/2.0) **2 # m^2

#Ndet_tng = np.array([1836, 950, 460])
Ndet_tng = np.array([1486, 638, 330])
fwhm_tng = np.array([25., 35., 50.])
A_beam_tng = 1.13 * fwhm_tng**2 # arcsec^2
crosspol = 1.0
det_eff = 1.0
map_eff = 1./np.sqrt(0.8)

# These values are for BLAST06, Laura's document
# They are the noise in a 5 hr, 1 deg^2 map for
# pol percentage error bars of 0.5 percent
noise_06 = np.array([1.44, 0.97, 0.52]) # (MJy/str) * sqrt(s)
noise_tng = noise_06 * np.sqrt((Ndet_06 * A_prim_06) /\
       (Ndet_tng * det_eff * A_prim_tng)) * crosspol\
        * map_eff # (MJy/str) * sqrt(s)

h = const.MKS_H
k = const.MKS_K
c = const.MKS_C
beta = 1.27

def arcsec2toDeg2(area):
    """ converts input map from arcsec^2 to deg^2 """
    return area / (3600**2)

def hour2sec(sec):
    """ converts hour to sec """
    return sec*3600

def obs_time(mapsize, sig_p, I_source, band_idx):
    """ Returns required observing time (hr) to make a map of map_size (deg^2)
        in band[band_idx], with an intensity in MJy/str of I_source,
        with pol fraction error bars of sig_p (percent) """
    n_beams = mapsize / arcsec2toDeg2(A_beam_tng[band_idx])
    #t_map = (2 * np.sqrt(2) * np.sqrt(n_beams) * noise_tng[band_idx] / 3600 * (I_source * sig_p))**2
    t_map = n_beams * (2.0 * np.sqrt(2.0) * (noise_tng[band_idx] / sig_p) / I_source)**2 # sec
    return t_map / 3600.

def I_min(mapsize, sig_p, t_map, band_idx):
    """ Returns minimum observable intensity (MJy/sr) achievable by BLAST-TNG
        by observing a map of map_size (square deg), with 1-sigma error bars in
        pol fraction, sig_p, in t_map hours."""
    n_beams = mapsize / arcsec2toDeg2(A_beam_tng[band_idx])
    Imin = 2.0 * np.sqrt(2.0) * (noise_tng[band_idx] / sig_p) *\
           np.sqrt(n_beams / hour2sec(t_map)) # MJy/str
    return Imin

def Bmod(nu, T):
    """ Returns the intensity from a source emitting as a modified BB with dust
        temperature T (K), and dust spectral index beta equal to 2."""
    one = nu**(beta)
    two = (2.0*h*nu**3 / c**2)
    three = 1./((np.exp((h*nu)/(k*T)) - 1.))
    return one*two*three

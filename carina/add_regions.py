import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import regions as reg

def car_regions(ax, f, wcs, ls = '-', add_text = True):
    regs = []

    # Southern bubble
    sb_coord = SkyCoord(160.50, -60.00, unit='deg', frame='fk5')
    #sb = reg.EllipseSkyRegion(center=sb_coord,
    #    height=30.0 * u.arcmin, width=15.0 * u.arcmin,
    #    angle=25 * u.deg).to_pixel(wcs)
    sb = reg.EllipseAnnulusSkyRegion(center=sb_coord,\
           inner_width=20 * u.arcmin,\
           outer_width=28 * u.arcmin,\
           inner_height=35 * u.arcmin,\
           outer_height=43 * u.arcmin,\
           angle= 25 * u.deg).to_pixel(wcs)
    regs.append(sb.as_artist(facecolor='none', edgecolor='white', linestyle = ls, lw=2))

    # TR 16
    ra_16 = Angle('10h45m10s').deg
    dec_16 = Angle('-59d43m0s').deg
    tr16_coord = SkyCoord(ra_16, dec_16, unit='deg', frame='fk5')
    tr16 = reg.CircleSkyRegion(center=tr16_coord, radius=5 * u.arcmin).to_pixel(wcs)
    regs.append(tr16.as_artist(facecolor='none', linestyle = ls, edgecolor='white', lw=2))

    # TR 14
    ra_14 = Angle('10h43m56s').deg
    dec_14 = Angle('-59d33m0s').deg
    tr14_coord = SkyCoord(ra_14, dec_14, unit='deg', frame='fk5')
    tr14 = reg.CircleSkyRegion(center=tr14_coord, radius=5 * u.arcmin).to_pixel(wcs)
    regs.append(tr14.as_artist(facecolor='none', edgecolor='white', linestyle = ls, lw=2))

    # Eta Car
    ra = Angle('10h45m3.546s').deg
    dec = Angle('-59d41m3.95s').deg
    etacar_coord = SkyCoord(ra, dec, unit='deg', frame='fk5')
    etacar = reg.CircleSkyRegion(center=etacar_coord, radius=1 * u.arcmin).to_pixel(wcs)
    regs.append(etacar.as_artist(facecolor='lime', edgecolor='lime', fill = True, lw=2))

    # Northern bubble
    nb_coord = SkyCoord(160.59, -59.23, unit='deg', frame='fk5')
    nb = reg.RectangleSkyRegion(center=nb_coord,
        height=30.0 * u.arcmin, width=50.0 * u.arcmin, angle = 1.2 * u.degree).to_pixel(wcs)
    regs.append(nb.as_artist(facecolor='none', linestyle = ls, edgecolor='white', lw=2))

    # Southern pillars
    sp_coord = SkyCoord(161.58, -59.98, unit='deg', frame='fk5')
    sp = reg.RectangleSkyRegion(center=sp_coord,
        height=20.0 * u.arcmin, width=40.0 * u.arcmin, angle = 0.6 *u.degree).to_pixel(wcs)
    regs.append(sp.as_artist(facecolor='none', linestyle = ls, edgecolor='white', lw=2))

    for r in regs:
        ax.add_patch(r)
    
    if add_text:
        f.add_label(161.58, -59.98, text = '5', size = 20, color = 'white')
        f.add_label(160.59, -59.23, text = '4', size = 20, color = 'white')
        f.add_label(160.51, -59.96, text = '1', size = 20, color = 'white')
        f.add_label(ra_14 - 0.1, dec_14, text = '2', size = 20, color = 'white')
        f.add_label(ra_16 + 0.08, dec_16, text = '3', size = 20, color = 'white')
    return

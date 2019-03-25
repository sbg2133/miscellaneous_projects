import numpy as np
from astropy.io import fits
import sys, os
from astropy.wcs import WCS
import aplpy
import glob
import matplotlib.pyplot as plt
import sens_calc as sens
save = 0
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ion()

rootdir = './galaxy_maps'
obj_list = ['m83', 'ngc1808', 'ngc4945', 'ngc2997', 'eso9713', 'ngc3621']
distance = np.array([8.9, 10.9, 3.9, 15.9, 4.2, 6.7])*1.0e6
scale = distance * np.tan(30./3600. * np.pi/180.)/1.0e3
print scale

def get_im_data(filepath):
    hdulist = fits.open(filepath)
    data = hdulist[1].data
    wcs = WCS(hdulist[1].header)
    return data, wcs

t_levels = np.linspace(1., 24., 10)
count = 0
for i in range(len(obj_list)):
    im = os.path.join(rootdir, obj_list[i], '250um.fits')
    data, wcs = get_im_data(im)
    cen_coord = wcs.wcs.crval
    t_req = sens.treq(0, 0.25, 0.05, data) # hours
    print "max t req =", np.nanmax(t_req)
    print "min t req =", np.nanmin(t_req)
    #fig = plt.figure(figsize = (15,12))
    fig = plt.figure()
    fig.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)
    #plt.imshow(t_req, vmin = 0.5, vmax = 24, cmap = 'inferno_r')
    hdu_t = fits.PrimaryHDU(data=t_req, header=wcs.to_header())
    ft = aplpy.FITSFigure(hdu_t, figure = fig, subplot=[0.1,0.1,0.45,0.8])
    ax = plt.gca()
    ax.set_facecolor('k')
    ft.recenter(cen_coord[0], cen_coord[1], width = 0.5, height = 0.5)
    if obj_list.index(obj_list[i]) == 0:
        vmin = 0.0
        vmax = 7.0
    if obj_list.index(obj_list[i]) == 1:
        vmin = 0.0
        vmax = 7.0
    if obj_list.index(obj_list[i]) == 2:
        vmin = 0.0
        vmax = 1.0
    if obj_list.index(obj_list[i]) == 3:
        vmin = 0.0
        vmax = 10.0
    if obj_list.index(obj_list[i]) == 4:
        vmin = 0.0
        vmax = 1.0
    if obj_list.index(obj_list[i]) == 5:
        vmin = 0.0
        vmax = 8.0
    ft.show_colorscale(vmin=vmin,vmax=vmax,stretch='linear', cmap = 'inferno_r')
    #levels = ft.show_contour(hdu_t, levels=10, filled = False, cmap = 'inferno_r', returnlevels = True) # hr
    ft.add_colorbar()
    ft.colorbar.set_location('top')
    ft.colorbar.set_pad(0.20)
    ft.colorbar.set_axis_label_font(size=16)
    ft.colorbar.set_font(size = 16)
    ft.colorbar.set_axis_label_text('t$_{req}$ [hr]')
    ft.tick_labels.set_yformat('dd.dd')
    ft.tick_labels.set_xformat('dd.dd')
    ft.axis_labels.set_font(size=16)
    ft.tick_labels.set_font(size=16)
    ft.ticks.set_color('white')
    ft.ticks.set_linewidth('2')
    #plt.imshow(t_req, vmin = 0.5, vmax = 24, cmap = 'inferno_r')
    hdu = fits.PrimaryHDU(data=data, header=wcs.to_header())
    f = aplpy.FITSFigure(hdu, figure = fig, subplot=[0.55,0.1,0.5,0.8])
    ax = plt.gca()
    ax.set_facecolor('k')
    f.recenter(cen_coord[0], cen_coord[1], width = 0.5, height = 0.5)
    f.show_grayscale()
    f.add_colorbar()
    f.colorbar.set_location('top')
    f.colorbar.set_pad(0.20)
    f.colorbar.set_axis_label_font(size=16)
    f.colorbar.set_font(size = 16)
    f.colorbar.set_axis_label_text(r'I$_{250}$ [MJy sr$^{-1}$]')
    #f.tick_labels.set_yformat('dd.dd')
    f.tick_labels.set_xformat('dd.dd')
    f.axis_labels.set_font(size=16)
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()
    f.tick_labels.set_font(size=16)
    f.ticks.set_color('white')
    f.ticks.set_linewidth('2')
    if save:
        plt.savefig(os.path.join(save_files_here, str(obj_list[count]) + '.eps'), format='eps', dpi = 500, bbox_inches = 'tight')
    count += 1

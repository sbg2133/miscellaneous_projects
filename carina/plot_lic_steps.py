import streamLines as sl
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
from getIQU import IQU
import scipy.ndimage as ndimage
from makePretty import pretty
import sys, os
from cycler import cycler
from matplotlib.collections import LineCollection
matplotlib.rcParams.update({'font.size':16})
save_files_here = "/home/wizwit/SESE_dissertation/figures/chapter6"
plt.ion()

band = sys.argv[1]
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']
pol_eff = [0.81, 0.79, 0.82]
blastpol_dir = './carinaData/smooth/5.0_arcmin_no_kernnorm'
filename = glob.glob(blastpol_dir + '/carinaneb_' + band + '_smoothed_5.0_rl_cal.fits')[0]
# load in I, Q, U for desired band
Ivals, Qvals, Uvals, pol_data, __ = IQU(filename, do_cov = True)
phi = np.deg2rad(pol_data[1][30:-30,260:-260])
p = pol_data[2][30:-30,260:-260]
p[p > 0.8] = np.nan
p[p < 0] = np.nan
sig_phi = np.deg2rad(pol_data[5][30:-30,260:-260])
sig_p = pol_data[5][30:-30,260:-260]
#phi2 = np.sqrt(phi**2 - sig_phi**2)
phi2 = phi - ((sig_phi**2 / 2.) / phi)
#I = Ivals[30:-30,260:-260][170:180,140:150]
#Q = Qvals[30:-30,260:-260][170:180,140:150]
#U = Uvals[30:-30,260:-260][170:180,140:150]
I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I
# Correct pvals as in Jamil's thesis, 5.7
pvals[pvals > 0.5] = np.nan
pvals[pvals < 0] = np.nan
pvals /= pol_eff[bands.index(band)]
phi = 0.5*np.arctan2(U,Q)
#dx = pvals*np.cos(phi)
#dy = pvals*np.sin(phi)
dx = np.cos(phi2)
dy = np.sin(phi2)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
xs, ys = np.meshgrid(X,Y)
xsize, ysize = len(X), len(Y)
vectors = np.array([dx,dy])
white = np.random.rand(xsize, ysize)

#  entire image, vectors downsampled x 5
fig = plt.figure()
ax = plt.gca()
q = sl.plot_vectors(ax, vectors, ys, xs, nskip = 5, alph = 1, col = 'k', pot = False)
#plt.axes().set_aspect('equal')
plt.xlabel("x [pix]")
plt.ylabel("y [pix]")
pretty()
plt.grid(False)
#plt.savefig(os.path.join(save_files_here, 'vectors_5.eps'), format='eps', dpi=1000, bbox_inches = 'tight')

#  entire image, vectors downsampled x 5
fig = plt.figure()
ax = plt.gca()
q = sl.plot_vectors(ax, vectors, ys, xs, nskip = 10, alph = 1, col = 'k', pot = False)
#plt.axes().set_aspect('equal')
plt.xlabel("x [pix]")
plt.ylabel("y [pix]")
pretty()
plt.grid(False)
#plt.savefig(os.path.join(save_files_here, 'vectors_10.eps'), format='eps', dpi=1000, bbox_inches = 'tight')
#  vectors and streamlines, entire image
fig = plt.figure()
ax = plt.gca()
c = ['firebrick']
steps = [21]
for i in range(len(steps)):
    print steps[i]
    sl.plot_streams(ax, vectors, xs, ys, nskip = 10, vec = False, Psteps = steps[i], col = c[i], alph = 0.7)
q = sl.plot_vectors(ax, vectors, ys, xs, nskip = 20, alph = 1, col = 'k', pot = False)
#hlines = np.column_stack(np.broadcast_arrays(X[0], Y, X[-1], Y))
#vlines = np.column_stack(np.broadcast_arrays(X, Y[0], X, Y[-1]))
#lines = np.concatenate([hlines, vlines]).reshape(-1, 2, 2)
#line_collection = LineCollection(lines, color="k", linewidths=1, alpha = 0.1)
#ax.add_collection(line_collection)
#plt.axes().set_aspect('equal')
#ax.set_xlim(140,150)
#ax.set_ylim(170, 180)
plt.xlabel("x [pix]")
plt.ylabel("y [pix]")
pretty()
#plt.pcolor(X, Y, I, alpha = 0.1, cmap='inferno')
plt.grid(False)
#fig.tight_layout()
#plt.savefig(os.path.join(save_files_here, 'vector_sl_21_ds10.eps'), format='eps', dpi=1000, bbox_inches = 'tight')

#  zoom
xlow = 265
xhigh = 280
ylow = 160
yhigh = 180
I = I[xlow:xhigh,ylow:yhigh]
Q = Q[xlow:xhigh,ylow:yhigh]
U = U[xlow:xhigh,ylow:yhigh]
phi3 = phi2[xlow:xhigh,ylow:yhigh]
pvals = pvals[xlow:xhigh,ylow:yhigh]
#Pvals = np.sqrt(Q**2 + U**2)
#pvals = Pvals/I
# Correct pvals as in Jamil's thesis, 5.7
pvals[pvals > 0.5] = np.nan
pvals[pvals < 0] = np.nan
pvals /= pol_eff[bands.index(band)]
#phi = 0.5*np.arctan2(U,Q)
#dx = pvals*np.cos(phi)
#dy = pvals*np.sin(phi)
dx = np.cos(phi3)
dy = np.sin(phi3)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
xs, ys = np.meshgrid(X,Y)
xsize, ysize = len(X), len(Y)
vectors = np.array([dx,dy])
white = np.random.rand(xsize, ysize)
fig = plt.figure()
ax = plt.gca()
c = ['firebrick']
steps = [3]
for i in range(len(steps)):
    print steps[i]
    sl.plot_streams(ax, vectors, xs, ys, nskip = 1, vec = False, Psteps = steps[i], col = c[i], alph = 0.9)
q = sl.plot_vectors(ax, vectors, ys, xs, nskip = 1, alph = 1, col = 'k', pot = False)
#hlines = np.column_stack(np.broadcast_arrays(X[0], Y, X[-1], Y))
#vlines = np.column_stack(np.broadcast_arrays(X, Y[0], X, Y[-1]))
#lines = np.concatenate([hlines, vlines]).reshape(-1, 2, 2)
#line_collection = LineCollection(lines, color="k", linewidths=1, alpha = 0.1)
#ax.add_collection(line_collection)
plt.axes().set_aspect('equal')
plt.xlabel("x [pix]")
plt.ylabel("y [pix]")
pretty()
#plt.pcolor(X, Y, I, alpha = 0.1, cmap='inferno')
plt.grid(False)
#fig.tight_layout()
#plt.savefig(os.path.join(save_files_here, 'vector_sl_zoom.eps'), format='eps', dpi=500, bbox_inches = 'tight')

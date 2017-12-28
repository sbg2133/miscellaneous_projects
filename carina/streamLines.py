import sys
import numpy as np
from getIQU import IQU
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage

#plt.ion()
#fig, ax = plt.subplots(figsize = (10.24, 7.68), dpi = 100)
#fig.set_facecolor("k")

def pix_idx(P, ys, xs):
    """Returns grid indices corresponding to point P"""
    x_idx = (np.abs(xs[0] - P[0])).argmin()
    y_idx = (np.abs(ys[:,0] - P[1])).argmin()
    return x_idx, y_idx

def get_vector(P, vectors, ys, xs):
    """Returns the vector located at pixel coordinate P"""
    x_idx, y_idx = pix_idx(P, ys, xs)
    vx, vy = vectors[0][y_idx, x_idx], vectors[1][y_idx, x_idx]
    angle = np.degrees(np.arctan2(vy, vx))
    return vx, vy, angle

def plot_pot():
    """plots the scalar field"""
    ax = plt.gca()
    im = ax.imshow(pot, origin = "lower", extent = [xs[0].min(), xs[0].max(),\
                ys[:,0].min(), ys[:,0].max()])
    return

def plot_vectors(vectors, ys, xs, nskip = 1, pot = False):
    """Creates an arrow plot of the vector field"""
    skip = (slice(None, None, nskip), slice(None, None, nskip))
    ax = plt.gca()
    if pot:
        plot_pot()
    mag = np.sqrt(vectors[0]**2 + vectors[1]**2)
    ax.quiver(xs[skip], ys[skip], (vectors[0]/mag)[skip], (vectors[1]/mag)[skip],\
             angles = 'xy', units = 'xy', scale_units = 'xy', headwidth = 2,\
             headaxislength = 2, headlength = 2, scale = 1)
    #ax.set(xticks = X, yticks = Y, aspect=1, title='Scratch', xlabel = 'x', ylabel = 'y')
    return

def new_P(idx, vectors, start, temp_pos, ys, xs, back = False):
    """Uses Euler's method to advect the streamline to the next position.
       @param idx: integers corresponding to start pixel position
       @param start: the starting coordinates (center of streamline)
       @param temp_pos: buffer for new advected position
       returns: temp_pos
                seg: distance to the previous temp_pos"""
    if (idx == 1):
        vx, vy = get_vector(start, vectors, ys, xs)[:2]
    else:
        vx, vy, angle = get_vector(temp_pos[idx - 1], vectors, ys, xs)
        if (np.isnan(vx) or np.isnan(vy)):
            pass
        #else:
            #vx2, vy2, angle2 = get_vector(temp_pos[idx - 2], vectors, ys, xs)
            #if (np.abs(angle2 - angle) >= 155.):
            #    return temp_pos, None
    if back:
        vx *= -1.
        vy *= -1.
    Px = np.int(temp_pos[idx - 1][0])
    Py = np.int(temp_pos[idx - 1][1])
    mag = np.sqrt(vx**2 + vy**2)
    s_top = ((Py + 1) - temp_pos[idx - 1][1])*(mag/vy)
    s_bot = (Py - temp_pos[idx - 1][0])*(mag/vy)
    s_right = ((Px + 1) - temp_pos[idx - 1][0])*(mag/vx)
    s_left = (Px - temp_pos[idx - 1][1])*(mag/vx)
    slist = np.array([s_top, s_bot, s_left, s_right])
    for s in slist:
        if (np.isnan(s)):
            return temp_pos, None
    if (slist < 0).all():
        s = np.min(np.abs(slist))
    else:
        s = np.min(slist[slist >= 0.])
    # add small amount to s to ensure that P is new pixel
    s += 0.08
    new_Px = temp_pos[idx - 1][0] + ((vx/mag)*s)
    new_Py = temp_pos[idx - 1][1] + ((vy/mag)*s)
    if (np.abs(new_Px - temp_pos[idx - 1][0]) > 2.):
        return temp_pos, None
    if (np.abs(new_Py - temp_pos[idx - 1][1]) > 2.):
        return temp_pos, None
    #if (new_Px > xs[0].max() or new_Px < xs[0].min()):
    #    return temp_pos, None
    #if (new_Py > ys[:,0].max() or new_Py < ys[:,0].min()):
    #    return temp_pos, None
    temp_pos.append((new_Px, new_Py))
    return temp_pos, s

def sl(start, vectors, ys, xs, plot = False):
    """Calculates a streamline centered on start"""
    # forward advection
    forward_pos = []
    forward_seg = []
    forward_pos.append(start)
    for i in range(1, Psteps):
        forward_pos, seg = new_P(i, vectors, start, forward_pos, ys, xs)
        if seg is not None:
            forward_seg.append(seg)
        else:
            break
    # backward advection
    back_pos = []
    back_seg = []
    back_pos.append(start)
    for i in range(1, Psteps):
        back_pos, seg = new_P(i, vectors, start, back_pos, ys, xs, back = True)
        if seg is not None:
            back_seg.append(seg)
        else:
            break
    streamline = list(reversed(forward_pos[1:]))
    streamline.extend(back_pos)
    # clean streamline of NaNs
    remove_idx = []
    count = 0
    for P in streamline:
        __, __, angle = get_vector(P, vectors, ys, xs)
        if (np.isnan(angle)):
            remove_idx.append(count)
        count += 1
    if remove_idx:
        for idx in sorted(remove_idx, reverse=True):
            del streamline[idx]
    #temp = clean_streamline(temp)
    streamline = np.array(streamline)
    if (plot):
        ax = plt.gca()
        for P in streamline:
            dx, dy, __ = get_vector(P, vectors, ys, xs)
            mag = np.sqrt(dx**2 + dy**2)
            ax.quiver(P[0], P[1], (dx/mag), (dy/mag), angles = 'xy',\
                units = 'xy', scale_units = 'xy', headwidth = 2, headaxislength = 2,\
                headlength = 2, scale = 1)
    if (plot):
        ax.scatter(streamline[:,0], streamline[:,1])
        ax.plot(streamline[:,0], streamline[:,1])
    return forward_seg, forward_pos, back_seg, back_pos, streamline

def plot_streams(ax, vectors, xs, ys, nskip = 1, vec = False, pot = False):
    """plots all streamlines. Launches a streamline from every
        grid point, modulo nskip.
        @param nskip: skip every nskip pixel
        @param vec: show arrow vectors
        @param pot: show potentials"""
    if vec:
        plot_vectors(vectors, ys, xs, nskip)
    if pot:
        plot_pot()
    for i in xs[0][::nskip]:
        for j in ys[:,0][::nskip]:
            __, __, __, __, streamline = sl([i, j], vectors, ys, xs)
            if not streamline.size:
                continue
            #if len(s.streamline[:,0]) < 5:
            #    continue
            ax.plot(streamline[:,0], streamline[:,1], color = 'white', alpha = 0.3)
    ax.set_facecolor("k")
    #ax.set(xticks = X, yticks = Y, aspect=1, title='Scratch', xlabel = 'x', ylabel = 'y')
    plt.tight_layout()
    return

def kern(k, s, L, hanning = False):
    """the convolution kernel"""
    # boxcar filter
    if not hanning:
        return k + s
    else:
        return k + (np.cos((s * np.pi) / L) + 1.)/2.

def partial_integral(forward_seg, forward_pos, back_seg, back_pos,\
                streamline, texture, ys, xs, hanning = False, back = False):
    """computes the line integral convolution in the forward
       or backward direction along a streamline.
       returns: Fsum - the LIC in one direction, for a single streamline
                hsum - a normalization factor"""
    if back:
        segs = back_seg
    else:
        segs = forward_seg
    np.insert(segs, 0, 0.)
    L = len(segs)
    klen = 31
    Fsum = 0.
    hsum = 0.
    s = 0.
    k0, k1 = 0., 0.
    for l in range(L):
        for i in range(1, len(segs) - 2):
            s += segs[i - 1]
            s_plus = s + segs[i + 1]
            if not hanning:
                k1 += s_plus
                k0 += s
            else:
                k1 += (np.cos((s_plus * np.pi) / klen) + 1.)/2.
                k0 += (np.cos((s * np.pi) / klen) + 1.)/2.
        h = k1 - k0
        hsum += h
        if back:
            tex_val = texture[pix_idx(back_pos[l], ys, xs)]
        else:
            tex_val = texture[pix_idx(forward_pos[l], ys, xs)]
        Fsum += tex_val * h
    return Fsum, hsum

def lic(forward_seg, forward_pos, back_seg, back_pos, streamline, texture, ys, xs):
    """performs a line integral convolution for a single streamline."""
    # compute forward integral
    # compute forward integral
    F_forward, h_forward = partial_integral(forward_seg, forward_pos,\
                back_seg, back_pos, streamline, texture, ys, xs, hanning = True)
    F_back, h_back = partial_integral(forward_seg, forward_pos,\
               back_seg, back_pos, streamline, texture, ys, xs, back = True, hanning = True)
    #print F_forward, h_forward
    #raw_input()
    if ((h_forward + h_back) == 0):
        temp_lic = 0.
    if ((F_forward + F_back) == 0):
        temp_lic = 0.
    else:
        temp_lic = (F_forward + F_back) / (h_forward + h_back)
    return temp_lic

def plot_lic(shape, vectors, texture):
    """plots the LIC"""
    xs, ys = np.meshgrid(X,Y)
    image = np.zeros((shape[0], shape[1]))
    ax = plt.gca()
    lics = []
    idx = 0
    for i in xs[0]:
        for j in ys[:,0]:
            start = [i,j]
            forward_seg, forward_pos, back_seg, back_pos,\
                                      streamline = sl(start, vectors, ys, xs)     
	    #if (start == [0,0]):
	        #print forward_pos
		#print 
		#print forward_seg
		#print 
		#print dx[0][0], dy[0][0]
            temp_lic = lic(forward_seg, forward_pos, back_seg, back_pos,\
            streamline, texture, ys, xs)
            lics.append(temp_lic)
            #print temp_lic
            #raw_input()
            if ((idx > 0) and temp_lic == 0.):
                temp_lic = lics[idx - 1]
            image[pix_idx(start, ys, xs)] = temp_lic
            idx += 1
    im = ax.imshow(image, origin="lower", cmap = "gist_heat")
    plt.autoscale(False)
    #ax.set(xticks = X, yticks = Y, aspect=1, title='Scratch', xlabel = 'x', ylabel = 'y')
    plt.tight_layout()
    return image

Psteps = 51
"""
band = sys.argv[1]
bands = ['250', '350', '500']
stokes = ['I', 'Q', 'U']
pol_eff = [0.81, 0.79, 0.82]
# load in I, Q, U for desired band
Ivals, Qvals, Uvals = IQU(band)
I = Ivals[30:-30,260:-260]
Q = Qvals[30:-30,260:-260]
U = Uvals[30:-30,260:-260]
Pvals = np.sqrt(Q**2 + U**2)
pvals = Pvals/I
# Correct pvals as in Jamil's thesis, 5.7
#pvals[pvals > 0.5] = np.nan
pvals /= pol_eff[bands.index(band)]
phi = 0.5*np.arctan2(U,Q)
dx = np.cos(phi)
dy = np.sin(phi)
mag = np.sqrt(dx**2 + dy**2)
X = np.linspace(0, I.shape[1], I.shape[1])
Y = np.linspace(0, I.shape[0], I.shape[0])
xs, ys = np.meshgrid(X,Y)
xsize, ysize = len(X), len(Y)
vectors = np.array([dx,dy])
white = np.random.rand(xsize, ysize)

#plot_streams(vectors, xs, ys, nskip = 10, vec = True)
#image = plot_lic([xsize, ysize], vectors, white)
"""

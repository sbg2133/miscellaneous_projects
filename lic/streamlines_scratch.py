import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from magnetic_dipole import dipole
plt.ion()
fig, ax = plt.subplots(figsize = (8, 5))

class Streamline:
    """A single streamline
       @param xstart: x,y starting position"""
    def __init__(self, xstart, ystart):
        self.start = (xstart, ystart)
        self.forward = []
        self.back = []
        self.s_forward = []
        self.s_back = []
        self.segs = []
        self.streamline = []

class VectorField():
    """A vector field"""
    #def __init__(self):

    def pix_idx(self, P):
        """Given (xcoord, ycoord), returns pixel indices"""
        x_idx = np.where(x[0] == np.int(P[0]))[0][0]
        y_idx = np.where(y[:,0] == np.int(P[1]))[0][0]
        return x_idx, y_idx

    def get_vector(self, P):
        """Returns the vector located at pixel coordinate P"""
        x_idx, y_idx = self.pix_idx(P)
        vx, vy = dx[y_idx, x_idx], dy[y_idx, x_idx]
        angle = np.degrees(np.arctan2(vy, vx))
        return vx, vy, angle

    def plot_pot(self):
        """plots the scalar field"""
        ax = plt.gca()
        im = ax.imshow(pot, origin = "lower", extent = [x[0].min(), x[0].max(),\
                    y[:,0].min(), y[:,0].max()])
        return

    def plot_vectors(self, nskip = 1, pot = False):
        """Creates an arrow plot of the vector field"""
        skip = (slice(None, None, nskip), slice(None, None, nskip))
        ax = plt.gca()
        if pot:
            plot_pot()
        mag = np.sqrt(dx**2 + dy**2)
        ax.quiver(x[skip], y[skip], (dx/mag)[skip], (dy/mag)[skip],\
                 angles = 'xy', units = 'xy', scale_units = 'xy', headwidth = 2,\
                 headaxislength = 2, headlength = 2, scale = 1)
        ax.set(xticks = (np.arange(x[0].min(), x[0].max(), 5)),\
              yticks = (np.arange(y[:,0].min(), y[:,0].max(), 5)),\
              aspect=1, title='Scratch', xlabel = 'x', ylabel = 'y')
        plt.margins(0,0)
        return

    def new_P(self, idx, start, temp_pos, back = False):
        """Uses Euler's method to advect the streamline to the next position.
           @param idx: integers corresponding to start pixel position
           @param start: the starting coordinates (center of streamline)
           @param temp_pos: buffer for new advected position
           returns: temp_pos
                    seg: distance to the previous temp_pos"""
        vector = self.get_vector
        if (idx == 1):
            vx, vy = vector(start)[:2]
        else:
            vx, vy, angle = vector(temp_pos[idx - 1])
            if (np.isnan(vx) or np.isnan(vy)):
                pass
            else:
                vx2, vy2, angle2 = vector(temp_pos[idx - 2])
                if (np.abs(angle2 - angle) >= 135.):
                    return temp_pos, None
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
        if (new_Px > x[0].max() or new_Px < x[0].min()):
            return temp_pos, None
        if (new_Py > y[:,0].max() or new_Py < y[:,0].min()):
            return temp_pos, None
        temp_pos.append((new_Px, new_Py))
        return temp_pos, s

    def streamline(self, xstart, ystart, plot = False):
        """a Streamline() instance
            @param xstart,ystart; x,y starting positions
            @param plot: If true, creates a plot of streamline
                with arrow vectors
            returns: a Streamline instance"""
        sl = Streamline(xstart, ystart)
        start = sl.start
        get_P = self.new_P
        vector = self.get_vector
        # forward advection
        forward_pos = []
        forward_seg = []
        segs = []
        forward_pos.append((xstart,ystart))
        for i in range(1, Psteps):
            forward_pos, seg = get_P(i, start, forward_pos)
            if seg is not None:
                forward_seg.append(seg)
                segs.append(seg)
            else:
                break
	sl.forward = forward_pos
	sl.s_forward = forward_seg
        # backward advection
        back_pos = []
        back_seg = []
        back_pos.append((xstart,ystart))
        for i in range(1, Psteps):
            back_pos, seg = get_P(i, start, back_pos, back = True)
            if seg is not None:
                back_seg.append(seg)
                segs.append(seg)
            else:
                break
        sl.back = back_pos
        sl.s_back = back_seg
        temp = list(reversed(forward_pos[1:]))
        temp.extend(back_pos)
        sl.segs = np.array(segs)
        # clean streamline of NaNs
        remove_idx = []
        count = 0
        for P in temp:
            __, __, angle = vector(P)
            if (np.isnan(angle)):
                remove_idx.append(count)
            count += 1
        if remove_idx:
            for idx in sorted(remove_idx, reverse=True):
                del temp[idx]
        #temp = self.clean_streamline(temp)
        sl.streamline = np.array(temp)
        if (plot):
            vector = self.get_vector
            ax = plt.gca()
            for P in sl.streamline:
                dx, dy, __ = self.get_vector(P)
                mag = np.sqrt(dx**2 + dy**2)
                ax.quiver(P[0], P[1], (dx/mag), (dy/mag), angles = 'xy',\
                    units = 'xy', scale_units = 'xy', headwidth = 2, headaxislength = 2,\
                    headlength = 2, scale = 1)
        if (plot):
            ax.scatter(sl.streamline[:,0], sl.streamline[:,1])
            ax.plot(sl.streamline[:,0], sl.streamline[:,1])
        return sl

    def plot_streams(self, nskip = 1, vectors = False, pot = False):
        """plots all streamlines. Launches a streamline from every
            grid point, modulo nskip.
            @param nskip: skip every nskip pixel
            @param vectors: f vectors, show arrow vectors
            @param pot: if pot, show potentials"""
        stream = self.streamline
        if vectors:
            self.plot_vectors()
        if pot:
            self.plot_pot()
        ax = plt.gca()
        for i in x[0][::nskip]:
            for j in y[:,0][::nskip]:
                s = stream(i, j)
                if not s.streamline.size:
                    continue
                #if len(s.streamline[:,0]) < 5:
                #    continue
                #ax.scatter(s.start[0], sl.start[1])
                ax.plot(s.streamline[:,0], s.streamline[:,1])
        return

    def kern(self, k, s, L, hanning = False):
        """the convolution kernel"""
        # boxcar filter
        if not hanning:
            return k + s
        else:
            return k + (np.cos((s * np.pi) / L) + 1.)/2.

    def partial_integral(self, sl, texture, hanning = False, back = False):
        """computes the line integral convolution in the forward
           or backward direction along a streamline.
           @param sl: a streamline instance
           @param texture: the background texture
           returns: Fsum - the LIC in one direction, for a single streamline
                    hsum - a normalization factor"""
        pix = self.pix_idx
        if back:
            segs = sl.s_back
        else:
            segs = sl.s_forward
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
                tex_val = texture[pix(sl.back[l])]
            else:
                tex_val = texture[pix(sl.forward[l])]
            Fsum += tex_val * h
        return Fsum, hsum

    def lic(self, sl, texture, nskip = 1):
        """performs a line integral convolution for a single streamline.
           @param sl: a streamline instance
           @param k: the convolution kernel"""
        # compute forward integral
        # compute forward integral
        F_forward, h_forward = self.partial_integral(sl, texture)
        F_back, h_back = self.partial_integral(sl, texture, back = True)
        if (np.sum(h_forward + h_back) == 0):
            temp_lic = 0.
        if ((F_forward + F_back) == 0):
            temp_lic = 0.
        else:
            temp_lic = (F_forward + F_back) / (h_forward + h_back)
        return temp_lic

    def plot_lic(self, texture, nskip = 1):
        """plots the LIC"""
        get_idx = self.pix_idx
        sl = self.streamline
        conv = self.lic
        self.image = np.zeros((pot.shape[0], pot.shape[1]))
        ax = plt.gca()
        ax.set_facecolor('black')
        lics = []
        idx = 0
        for i in x[0][::nskip]:
            for j in y[:,0][::nskip]:
                s = sl(i, j)
                temp_lic = conv(s, texture, nskip = nskip)
                lics.append(temp_lic)
                if ((idx > 0) and temp_lic == 0.):
                    temp_lic = lics[idx - 1]
                self.image[get_idx(s.start)] = temp_lic
                idx += 1
        im = ax.imshow(self.image, origin="lower", cmap = "gist_heat")
        plt.tight_layout()
        return

xsize, ysize = 25, 25
y, x = np.mgrid[0:2*ysize, 0:2*xsize]
Psteps = 25

"""
### point masses ###
pot = np.zeros((2*xsize, 2*ysize))
mass = [2**10, 2**10]
pos = [(15.4,15.2), (36.8,39.1)]
for i in range(len(pos)):
   r = np.sqrt((x - pos[i][0])**2 + (y - pos[i][1])**2)
   pot += mass[i] / r
pot[~np.isfinite(pot)] = 0.0
interp = 4
if (interp > 1):
    pot = ndimage.zoom(pot, interp, order = 1)
    y, x = np.mgrid[:pot.shape[0], :pot.shape[1]]
    dy, dx = np.gradient(pot, np.diff(y[:2,0])[0], np.diff(x[0,:2])[0])
"""
### magnetic dipole ###
Bx, By = dipole(m=[0.3, 0.1], r=np.mgrid[0:2*ysize, 0:2*xsize], r0=[25.1,25.1])
pot = np.zeros((2*xsize, 2*ysize))
interp = 4
if (interp > 1):
    rx = (interp * 100 / 4.) + 0.1
    ry = (interp * 100 / 4.) + 0.1
    pot = ndimage.zoom(pot, interp, order = 1)
    y, x = np.mgrid[:pot.shape[0], :pot.shape[1]]
    dx, dy = dipole(m=[0.3, 0.1], r=np.mgrid[0:pot.shape[0], 0:pot.shape[1]],\
                         r0=[rx,ry])
white = np.random.rand(pot.shape[0], pot.shape[1])
scratch = VectorField()
#scratch.plot_lic(white)


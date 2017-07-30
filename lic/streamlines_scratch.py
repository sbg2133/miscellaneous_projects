import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
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
        self.angles = []

class VectorField:
    """A vector field
       @param xsize, ysize: x,y dimensions
       @param potentials: 2-D array of scalar potentials
       @param dx, dy: gradient of potentials in x,y direction"""
    def __init__(self, x, y, dx, dy, potentials, interp = 1):
        self.pot = potentials
        self.y, self.x = y, x
        self.dx, self.dy = dx, dy
        del dx
        del dy
        del x
        del y
        del potentials
        self.Psteps = 30
        if (interp > 1):
            self.pot = ndimage.zoom(self.pot, interp, order = 1)
            self.y, self.x = np.mgrid[:self.pot.shape[0], :self.pot.shape[1]]
            self.white = np.random.rand(self.pot.shape[0], self.pot.shape[1])
            self.dy, self.dx = np.gradient(self.pot, np.diff(self.y[:2,0])[0],\
                      np.diff(self.x[0,:2])[0])
        self.white = np.random.rand(self.pot.shape[0], self.pot.shape[1])

    def pix_idx(self, P):
        """Given (xcoord, ycoord), returns pixel indices"""
        x_idx = np.where(self.x[0] == np.int(P[0]))[0][0]
        y_idx = np.where(self.y[:,0] == np.int(P[1]))[0][0]
        return x_idx, y_idx

    def get_vector(self, P):
        """Given pixel coord (x, y), returns vector (dx, dy)"""
        x_idx, y_idx = self.pix_idx(P)
        dx, dy = self.dx[y_idx, x_idx], self.dy[y_idx, x_idx]
        angle = np.degrees(np.arctan2(dy, dx))
        return dx, dy, angle

    def plot_pot(self):
        """plots the scalar field"""
        ax = plt.gca()
        im = ax.imshow(pot, origin = "lower", extent = [self.x[0].min(),\
                 self.x[0].max(), self.y[:,0].min(),\
                 self.y[:,0].max()])
        return

    def plot_vectors(self, nskip = 1, pot = False):
        """Creates an arrow plot of the vector field"""
        skip = (slice(None, None, nskip), slice(None, None, nskip))
        ax = plt.gca()
        if pot:
            plot_pot()
        mag = np.sqrt(self.dx**2 + self.dy**2)
        ax.quiver(self.x[skip], self.y[skip], (self.dx/mag)[skip], (self.dy/mag)[skip],\
                 angles = 'xy', units = 'xy', scale_units = 'xy', headwidth = 2,\
                 headaxislength = 2, headlength = 2, scale = 1)
        ax.set(xticks = (np.arange(self.x[0].min(), self.x[0].max(), 5)),\
              yticks = (np.arange(self.y[:,0].min(), self.y[:,0].max(), 5)),\
              aspect=1, title='Scratch', xlabel = 'x', ylabel = 'y')
        plt.margins(0,0)
        return

    def new_P(self, idx, start, temp_pos, back = False):
        """Given a pixel index int(x,y) and streamline instance (sl),
           uses Euler's method to advect the streamline to the next
           position.
           @param idx: integers corresponding to start pixel position
           @param temp_stream: a list to hold the advected positions
           returns: -1 on failure, 0 on success"""
        vector = self.get_vector
        if (idx == 1):
            dx, dy = vector(start)[:2]
        else:
            dx, dy, angle = vector(temp_pos[idx - 1])
            if (np.isnan(dx) or np.isnan(dy)):
                pass
            else:
                dx2, dy2, angle2 = vector(temp_pos[idx - 2])
                if (np.abs(angle2 - angle) >= 135.):
                    return temp_pos, None
        if back:
            dx *= -1.
            dy *= -1.
        x = np.int(temp_pos[idx - 1][0])
        y = np.int(temp_pos[idx - 1][1])
        mag = np.sqrt(dx**2 + dy**2)
        s_top = ((y + 1) - temp_pos[idx - 1][1])*(mag/dy)
        s_bot = (y - temp_pos[idx - 1][0])*(mag/dy)
        s_right = ((x + 1) - temp_pos[idx - 1][0])*(mag/dx)
        s_left = (x - temp_pos[idx - 1][1])*(mag/dx)
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
        new_Px = temp_pos[idx - 1][0] + ((dx/mag)*s)
        new_Py = temp_pos[idx - 1][1] + ((dy/mag)*s)
        if (np.abs(new_Px - temp_pos[idx - 1][0]) > 2.):
            return temp_pos, None
        if (np.abs(new_Py - temp_pos[idx - 1][1]) > 2.):
            return temp_pos, None
        if (new_Px > self.x[0].max() or new_Px < self.x[0].min()):
            return temp_pos, None
        if (new_Py > self.y[:,0].max() or new_Py < self.y[:,0].min()):
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
        # forward advection
        forward_pos = []
        forward_seg = []
        segs = []
        forward_pos.append((xstart,ystart))
        for i in range(1, self.Psteps):
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
        for i in range(1, self.Psteps):
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
        temp = self.clean_streamline(temp)
        sl.streamline = temp
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

    def clean_streamline(self, sl_temp):
        """removes NaN values from streamline
            @param sl: a list of streamline values"""
        remove_idx = []
        count = 0
        for P in sl_temp:
            dx, dy, angle = self.get_vector(P)
            if (np.isnan(angle)):
                remove_idx.append(count)
            count += 1
        if remove_idx:
            for idx in sorted(remove_idx, reverse=True):
                del sl_temp[idx]
        return np.array(sl_temp)

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
        for x in self.x[0][::nskip]:
            for y in self.y[:,0][::nskip]:
                s = streamline(x, y)
                if not s.streamline.size:
                    continue
                #if len(s.streamline[:,0]) < 5:
                #    continue
                #ax.scatter(s.start[0], sl.start[1])
                ax.plot(s.streamline[:,0], s.streamline[:,1])
        return

    def kern(self, k, s, L, hanning = False):
        # boxcar filter
        if not hanning:
            return k + s
        else:
            return k + (np.cos((s * np.pi) / L) + 1.)/2.

    def partial_integral(self, sl, texture, back = False):
        kern = self.kern
        pix = self.pix_idx
        if back:
            segs = sl.s_back
        else:
            segs = sl.s_forward
        np.insert(segs, 0, 0.)
        L = len(segs)
        Fsum = 0.
        hsum = 0.
        s = 0.
        k0, k1 = 0., 0.
        for l in range(L):
            for i in range(1, len(segs) - 2):
                s += segs[i - 1]
                s_plus = s + segs[i + 1]
                k1 = kern(k1, s_plus, 256, hanning = True)
                k0 = kern(k0, s, 256, hanning = True)
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
        get_idx = self.pix_idx
        sl = self.streamline
        conv = self.lic
        self.image = np.zeros((self.pot.shape[0], self.pot.shape[1]))
        ax = plt.gca()
        ax.set_facecolor('black')
        lics = []
        idx = 0
        for x in self.x[0][::nskip]:
            for y in self.y[:,0][::nskip]:
                s = sl(x, y)
                temp_lic = conv(s, texture, nskip = nskip)
                lics.append(temp_lic)
                if ((idx > 0) and temp_lic == 0.):
                    temp_lic = lics[idx - 1]
                self.image[get_idx(s.start)] = temp_lic
                idx += 1
        im = ax.imshow(self.image, origin="lower", cmap = "gist_heat")
                # extent = [self.x[0].min(),\
                # self.x[0].max(), self.y[:,0].min(),\
                # self.y[:,0].max()])
        return

xsize, ysize = 25, 25
y, x = np.mgrid[0:2*ysize, 0:2*xsize]
#G = 6.67e-11 # m^3 / s^2 * kg
#m_earth = 6.0e24 # Earth mass, kg
# Number of point masses
pot = np.zeros((2*xsize, 2*ysize))
mass = [2**10, 2**10]
pos = [(15.4,15.2), (36.8,39.1)]
for i in range(len(pos)):
   r = np.sqrt((x - pos[i][0])**2 + (y - pos[i][1])**2)
   #r[r == 0.0] = np.nan
   pot += mass[i] / r
pot[~np.isfinite(pot)] = 0.0
dy, dx = np.gradient(pot, np.diff(y[:3,0])[0], np.diff(x[0,:2])[0])
scratch = VectorField(x, y, dx, dy, pot, interp = 4)
scratch.plot_lic(scratch.white)


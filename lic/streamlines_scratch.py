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
        self.vectors = []
        self.angles = []
        self.lic = []

class VectorField:
    """A vector field
       @param xsize, ysize: x,y dimensions
       @param potentials: 2-D array of scalar potentials
       @param dx, dy: gradient of potentials in x,y direction"""
    def __init__(self, xsize, ysize, dx, dy, potentials, interp = 1):
        self.pot = potentials
        self.white = np.random.rand(self.pot.shape[0], self.pot.shape[1])
        self.y, self.x = np.mgrid[0:2*ysize, 0:2*xsize]
        self.dx, self.dy = dx, dy
        self.Psteps = 30
        if (interp > 1):
            self.pot = ndimage.zoom(self.pot, interp, order = 1)
            self.y, self.x = np.mgrid[:self.pot.shape[0], :self.pot.shape[1]]
            self.white = np.random.rand(self.pot.shape[0], self.pot.shape[1])
            self.dy, self.dx = np.gradient(self.pot, np.diff(self.y[:2,0])[0],\
                      np.diff(self.x[0,:2])[0])
    def pix_idx(self, P):
        """Given (xcoord, ycoord), returns pixel indices"""
        x_idx = np.where(self.x[0] == np.int(P[0]))[0][0]
        y_idx = np.where(self.y[:,0] == np.int(P[1]))[0][0]
        return x_idx, y_idx

    def get_vector(self, P, norm = True):
        """Given pixel coord (x, y), returns vector (dx, dy)"""
        x_idx, y_idx = self.pix_idx(P)
        dx, dy = self.dx[y_idx, x_idx], self.dy[y_idx, x_idx]
        angle = np.degrees(np.arctan2(dy, dx))
        mag = np.sqrt(dx**2 + dy**2)
        if (norm == False):
            return dx, dy, angle
        return dx/mag, dy/mag, angle

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
        for i in range(self.pot.shape[0]):
            for j in range(self.pot.shape[1]):
                P = (self.x[i][j], self.y[i][j])
                dx, dy, __ = self.get_vector(P)
                ax.arrow(self.x[i][j], self.y[i][j], dx, dy,\
                           head_width=0.05, head_length=0.1, fc='k', ec='k')
        ax.set(xticks = (np.arange(self.x[0].min() - 1, self.x[0].max() + 2, 5)),\
              yticks = (np.arange(self.y[:,0].min() - 1, self.y[:,0].max() + 2, 5)),\
              aspect=1, title='Scratch', xlabel = 'x', ylabel = 'y')
        plt.tight_layout()
        return

    def new_P(self, idx, sl, back = False):
        """Given a pixel index int(x,y) and streamline instance (sl),
           uses Euler's method to advect the streamline to the next
           position.
           @param idx: integers corresponding to start pixel position
           @param sl: a Streamline() instance
           returns: -1 on failure, 0 on success"""
        if back:
            line = sl.back
        else:
            line = sl.forward
        if (idx == 1):
            dx, dy = self.get_vector(sl.start, norm = False)[:2]
        else:
            dx, dy, angle = self.get_vector(line[idx - 1], norm = False)
            if (np.isnan(dx) or np.isnan(dy)):
                #line.append((np.nan, np.nan))
                pass
            else:
                dx2, dy2, angle2 = self.get_vector(line[idx - 2], norm = False)
                if (np.abs(angle2 - angle) >= 135.):
                    return -1
        if back:
            dx *= -1.
            dy *= -1.
        x = np.int(line[idx - 1][0])
        y = np.int(line[idx - 1][1])
        mag = np.sqrt(dx**2 + dy**2)
        s_top = ((y + 1) - line[idx - 1][1])*(mag/dy)
        s_bot = (y - line[idx - 1][0])*(mag/dy)
        s_right = ((x + 1) - line[idx - 1][0])*(mag/dx)
        s_left = (x - line[idx - 1][1])*(mag/dx)
        slist = np.array([s_top, s_bot, s_left, s_right])
        for s in slist:
            if (np.isnan(s)):
                return -1
        if (slist < 0).all():
            s = np.min(np.abs(slist))
        else:
            s = np.min(slist[slist >= 0.])
        # add small amount to s to ensure that P is new pixel
        s += 0.08
        new_Px = line[idx - 1][0] + ((dx/mag)*s)
        new_Py = line[idx - 1][1] + ((dy/mag)*s)
        if (np.abs(new_Px - line[idx - 1][0]) > 2.):
            return -1
        if (np.abs(new_Py - line[idx - 1][1]) > 2.):
            return -1
        if (new_Px > self.x[0].max() or new_Px < self.x[0].min()):
            return -1
        if (new_Py > self.y[:,0].max() or new_Py < self.y[:,0].min()):
            return -1
        if back:
            sl.back.append((new_Px, new_Py))
            sl.s_back.append(s)
        else:
            sl.forward.append((new_Px, new_Py))
            sl.s_forward.append(s)
        sl.segs.append(s)
        return 0

    def streamline(self, xstart, ystart, plot = False):
        """a Streamline() instance
            @param xstart,ystart; x,y starting positions
            @param plot: If true, creates a plot of streamline
                with arrow vectors
            returns: a Streamline instance"""
        sl = Streamline(xstart, ystart)
        sl.forward.append((xstart, ystart))
        sl.back.append((xstart, ystart))
        # forward advection
        for i in range(1, self.Psteps):
            if (self.new_P(i, sl) < 0):
                break
        # backward advection
        for i in range(1, self.Psteps):
            if (self.new_P(i, sl, back = True) < 0):
                break
        temp = list(reversed(sl.forward[1:]))
        for (x,y) in sl.back:
            temp.append((x,y))
        sl.streamline = temp
        for P in sl.streamline:
            dx, dy, angle = self.get_vector(P)
            sl.vectors.append((dx,dy))
            sl.angles.append(angle)
        sl.vectors = np.array(sl.vectors)
        sl.angles = np.array(sl.angles)
        sl.segs = np.array(sl.segs)
        self.clean_streamline(sl)
        if (plot):
            ax = plt.gca()
            for P in sl.streamline:
                ax.arrow(P[0], P[1], dx, dy,\
                head_width=0.05, head_length=0.1, fc='k', ec='k')
        if (plot):
            ax.scatter(sl.streamline[:,0], sl.streamline[:,1])
            ax.plot(sl.streamline[:,0], sl.streamline[:,1])
        return sl

    def clean_streamline(self, sl):
        """removes NaN values from streamline
            @param sl: a streamline instance"""
        remove_idx = []
        for i in range(len(sl.streamline)):
            if (np.isnan(sl.angles[i])):
                remove_idx.append(i)
        sl.angles = sl.angles[~np.isnan(sl.angles)]
        sl.vectors = sl.vectors[~np.isnan(sl.vectors)]
        if remove_idx:
            for idx in sorted(remove_idx, reverse=True):
                del sl.streamline[idx]
        sl.streamline = np.array(sl.streamline)
        return

    def plot_streams(self, nskip = 1, vectors = False, pot = False):
        """plots all streamlines. Launches a streamline from every
            grid point, modulo nskip.
            @param nskip: skip every nskip pixel
            @param vectors: f vectors, show arrow vectors
            @param pot: if pot, show potentials"""
        if vectors:
            self.plot_vectors()
        if pot:
            self.plot_pot()
        ax = plt.gca()
        for x in self.x[0][::nskip]:
            for y in self.y[:,0][::nskip]:
                s = self.streamline(x, y)
                if not s.streamline.size:
                    continue
                if len(s.streamline[:,0]) < 5:
                    continue
                #ax.scatter(s.start[0], sl.start[1])
                ax.plot(s.streamline[:,0], s.streamline[:,1])
        return

    def partial_integral(self, sl, kernel, texture, back = False):
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
                k1 += s_plus
                k0 += s
            h = k1 - k0
            hsum += h
            if back:
                tex_val = texture[self.pix_idx(sl.back[l])]
            else:
                tex_val = texture[self.pix_idx(sl.forward[l])]
            Fsum += tex_val * h
        return Fsum, hsum

    def lic(self, sl, texture, nskip = 1):
        """performs a line integral convolution for a single streamline.
           @param sl: a streamline instance
           @param k: the convolution kernel"""
        # Boxcar kernel
        klen = 31
        kernel = np.ones(klen)
        #k = np.sin(np.arange(klen)*np.pi/klen
        # compute forward integral
        # compute forward integral
        F_forward, h_forward = self.partial_integral(sl, kernel, texture)
        F_back, h_back = self.partial_integral(sl, kernel, texture, back = True)
        if (np.sum(h_forward + h_back) == 0):
            sl.lic = 0.
            #sl.lic = texture[self.pix_idx(sl.start)]
        if ((F_forward + F_back) == 0):
            sl.lic = 0.
            #sl.lic = texture[self.pix_idx(sl.start)]
        else:
            sl.lic = (F_forward + F_back) / (h_forward + h_back)
        #self.lics.append(sl.lic)
        return

    def plot_lic(self, texture, nskip = 1):
        self.image = np.zeros((self.pot.shape[0], self.pot.shape[1]))
        ax = plt.gca()
        ax.set_facecolor('black')
        lics = []
        idx = 0
        for x in self.x[0][::nskip]:
            for y in self.y[:,0][::nskip]:
                s = self.streamline(x, y)
                self.lic(s, texture, nskip = nskip)
                lics.append(s.lic)
                if ((idx > 0) and s.lic == 0.):
                    s.lic = lics[idx - 1]
                self.image[self.pix_idx(s.start)] = s.lic
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
mass = [2**5, 2**5, 2**5, 2**5, 2**5, 2**5]
pos = [(15,15), (36,96), (20,70), (40,40), (73, 53)]
#pos = [[30., 30.]]
#mass = [[2**20]]
for i in range(len(pos)):
   r = np.sqrt((x - pos[i][0])**2 + (y - pos[i][1])**2)
   r[r == 0.0] = np.nan
   pot += mass[i] / r
dy, dx = np.gradient(pot, np.diff(y[:2,0])[0], np.diff(x[0,:2])[0])
scratch = VectorField(pot.shape[0]/2, pot.shape[0]/2, dx, dy, pot, interp = 8)
#scratch.plot_lic(white)


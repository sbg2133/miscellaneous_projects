from subprocess import call
import sys
import numpy as np
import matplotlib.pyplot as plt
from magnetic_dipole import dipole
plt.ion()
plt.figure(figsize = (10.24, 7.68), dpi = 100)

xsize, ysize = int(sys.argv[1]), int(sys.argv[2])
xmax, ymax = 200, 200
X = np.linspace(0, xmax, xsize)
Y = np.linspace(0, ymax, ysize)
x, y = np.meshgrid(X,Y)

### magnetic dipole ###
dx, dy = dipole(m=[5., 5.], r=np.meshgrid(X,Y), r0=[xmax/2. + 0.1, ymax/2. + 0.3]).astype('float32')
vectors = np.array([dx,dy])
white = np.random.rand(xsize, ysize)
with file('texture.dat', 'w') as outfile:
    for row in white:
        np.savetxt(outfile, row, newline = " ")
	outfile.write('\n')
with file('dx.dat', 'w') as outfile:
    for row in dx:
        np.savetxt(outfile, row, newline = " ")
	outfile.write('\n')
with file('dy.dat', 'w') as outfile:
    for row in dy:
        np.savetxt(outfile, row, newline = " ")
	outfile.write('\n')

command = ["./lic", str(xsize), str(ysize)]
call(command)

lic = np.loadtxt("./lic.dat")
plt.imshow(lic, cmap = "viridis", interpolation = "sinc")
plt.tight_layout()


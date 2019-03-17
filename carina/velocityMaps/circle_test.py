import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 32)
y = np.arange(0, 32)
arr = np.zeros((y.size, x.size))

cx = 12.
cy = 16.
r = 5.

# The two lines below could be merged, but I stored the mask
# for code clarity.
mask = (x[np.newaxis,:]-cx)**2 + (y[:,np.newaxis]-cy)**2 < r**2
arr[mask] = 123.

# This plot shows that only within the circle the value is set to 123.
plt.figure()
plt.pcolormesh(x, y, arr)
plt.colorbar()
plt.show()

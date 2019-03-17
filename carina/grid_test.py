import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

X = np.linspace(0, 410, 410)
Y = np.linspace(0, 410, 410)

from matplotlib.collections import LineCollection


plt.figure(figsize=(12, 7))

hlines = np.column_stack(np.broadcast_arrays(X[0], Y, X[-1], Y))
vlines = np.column_stack(np.broadcast_arrays(X, Y[0], X, Y[-1]))
lines = np.concatenate([hlines, vlines]).reshape(-1, 2, 2)
line_collection = LineCollection(lines, color="k", linewidths=1, alpha = 0.5)
ax = plt.gca()
ax.add_collection(line_collection)
ax.set_xlim(X[0], X[-1])
ax.set_ylim(Y[0], Y[-1])
plt.show()

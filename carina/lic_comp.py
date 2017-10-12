import numpy as np
import matplotlib.pyplot as plt

plt.ion()
lic250 = np.load('./lic250.npy')
lic350 = np.load('./lic350.npy')
lic500 = np.load('./lic500.npy')

f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, dpi = 100)
ax1.imshow(lic250, cmap = "inferno", interpolation = "gaussian")
ax2.imshow(lic350, cmap = "Greens", interpolation = "gaussian")
ax3.imshow(lic500, cmap = "Blues", interpolation = "gaussian")
plt.tight_layout()


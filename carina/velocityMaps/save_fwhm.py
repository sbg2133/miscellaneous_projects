import numpy as np
import sys
from fitFWHM import vfwhm

fit_type = sys.argv[1]
data = np.load('./large_files/v_array_skip3.npy')

# Lorentzian
if fit_type == 'lor':
    fwhm = vfwhm(data)
    np.save('./large_files/fwhm_HI_lor.npy', fwhm)
# Gaussian
else:
    fwhm = vfwhm(data, lor = False)
    np.save('./large_files/fwhm_HI_gauss.npy', fwhm)



import numpy as np
import matplotlib.pyplot as plt

"""Temporarily borrowed from:
   http://michal.rawlik.pl/2015/03/12/magnetic-dipole-in-python/
   Thanks!"""

def dipole(m, r, r0):
    R = np.subtract(np.transpose(r), r0).T
    norm_R = np.sqrt(np.einsum("i...,i...", R, R))
    m_dot_R = np.tensordot(m, R, axes=1)
    B = 3 * m_dot_R * R / norm_R**5 - np.tensordot(m, 1 / norm_R**3, axes=0)
    #B *= 1e-7
    return B

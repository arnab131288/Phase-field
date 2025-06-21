"""
Compute anisotropy factors and mixed derivative for a point.
"""
import math
from numba import njit
from config import DELTA, FOLD, rdX, rdY

@njit
def compute_anisotropy(phi, i, j):
    dpx = (phi[i+1, j] - phi[i-1, j]) * 0.5 * rdX
    dpy = (phi[i, j+1] - phi[i, j-1]) * 0.5 * rdY
    dxy = (
        phi[i+1, j+1] - phi[i+1, j-1]
        - phi[i-1, j+1] + phi[i-1, j-1]
    ) * 0.25 * rdX * rdY
    theta = math.atan2(dpy, dpx)
    eta  = 1 + DELTA * math.cos(FOLD * theta)
    eta1 = -FOLD * DELTA * math.sin(FOLD * theta)
    eta2 = -FOLD * FOLD * DELTA * math.cos(FOLD * theta)
    return eta, eta1, eta2, theta, dxy




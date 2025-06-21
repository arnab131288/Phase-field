"""
Compute Laplacian stencil at a grid point using finite differences.
"""
import math
from numba import njit
from config import dx2, dy2

@njit
def compute_laplacians(f, i, j):
    lapX = (f[i+1, j] - 2*f[i, j] + f[i-1, j]) * dx2
    lapY = (f[i, j+1] - 2*f[i, j] + f[i, j-1]) * dy2
    return lapX, lapY

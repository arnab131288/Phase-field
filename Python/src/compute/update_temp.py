"""
Update temperature field with latent heat source from phase change.
"""
from numba import njit
from config import DELTA_T, K

@njit
def update_temp(temp_old, lapX, lapY, phi_new, phi_old):
    lap = lapX + lapY
    source = (phi_new - phi_old) / DELTA_T * K
    return temp_old + (lap + source) * DELTA_T



"""
Compute phase-field driving force with thermal noise.
"""
import math
import random
from numba import njit
from config import A, ALPHA, GAMMA, T_E

@njit
def compute_driving(phi_ij, temp_ij):
    noise = A * (random.random() - 0.5)
    m = (ALPHA / math.pi) * math.atan(GAMMA * (T_E - temp_ij))
    return phi_ij * (1 - phi_ij) * (phi_ij - 0.5 + m + noise)



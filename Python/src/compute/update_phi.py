"""
Update phase field at a grid point including anisotropic effects.
"""
import math
from numba import njit
from config import EPSILON, DELTA_T, TAU

@njit
def update_phi(phi_old, lapX, lapY, driving, eta, eta1, eta2, theta, dxy):
    lap = lapX + lapY
    c1 = lapY - lapX
    s2t = math.sin(2 * theta)
    c2t = math.cos(2 * theta)
    term_aniso = (EPSILON**2 * eta * eta1 * (s2t * c1 + 2 * c2t * dxy)
                  - 0.5 * EPSILON**2 * (eta1**2 + eta * eta2)
                  * (2 * s2t * dxy - lap - c2t * c1))
    total = EPSILON**2 * lap + driving + term_aniso
    return phi_old + total * (DELTA_T / TAU)




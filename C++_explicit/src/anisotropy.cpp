#include "header.hpp"
#include <cmath>

// Macro to compute flattened array index for 3D data
#define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Compute anisotropy functions and their derivatives for the phase-field.
 *
 * For each interior grid point, computes:
 *   theta = atan2(DERY_c, DERX_c)
 *   a_c = epsilon* [1 + delta cos(j (theta - theta0))]
 *   a_p = -epsilon * delta * j * sin(j (theta - theta0))
 * Then computes the same for neighboring derivative directions.
 *
 * @param fb      FieldBuffers containing gradient arrays and outputs for anisotropy
 * @param params  Simulation parameters including epsilon, delta, j, theta0, grid dims
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1]
 */
void computeAnisotropy(FieldBuffers *fb, const SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    const double eps   = params->epsilon;
    const double delta = params->delta;
    const int    jmult = params->j;
    const double theta0= params->theta_0;

    // Loop over interior grid (assuming k=0 for 2D or first layer for 3D)
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            int idx = IDX(i, j, (dim==3)?1:0);

            // Compute interface normal angle
            double theta = std::atan2(fb->DERY_c[idx], fb->DERX_c[idx]);
            // Anisotropy and its derivative
            fb->ac[idx]   = eps * (1.0 + delta * std::cos(jmult * (theta - theta0)));
            fb->ac_p[idx] = -eps * (delta * jmult * std::sin(jmult * (theta - theta0)));

            // Neighbors' interface angles
            double theta_r = std::atan2(fb->DERY_right[idx],  fb->DERX_right[idx]);
            double theta_l = std::atan2(fb->DERY_left[idx],   fb->DERX_left[idx]);
            double theta_t = std::atan2(fb->DERY_top[idx],    fb->DERX_top[idx]);
            double theta_b = std::atan2(fb->DERY_bottom[idx], fb->DERX_bottom[idx]);

            // Anisotropy at neighbors
            fb->ac_right[idx]  = eps * (1.0 + delta * std::cos(jmult * (theta_r - theta0)));
            fb->ac_left[idx]   = eps * (1.0 + delta * std::cos(jmult * (theta_l - theta0)));
            fb->ac_top[idx]    = eps * (1.0 + delta * std::cos(jmult * (theta_t - theta0)));
            fb->ac_bottom[idx] = eps * (1.0 + delta * std::cos(jmult * (theta_b - theta0)));

            // Derivative of anisotropy at neighbors
            fb->ac_p_right[idx]  = -eps * (delta * jmult * std::sin(jmult * (theta_r - theta0)));
            fb->ac_p_left[idx]   = -eps * (delta * jmult * std::sin(jmult * (theta_l - theta0)));
            fb->ac_p_top[idx]    = -eps * (delta * jmult * std::sin(jmult * (theta_t - theta0)));
            fb->ac_p_bottom[idx] = -eps * (delta * jmult * std::sin(jmult * (theta_b - theta0)));
        }
    }
}

#undef IDX






#include "header.hpp"

// Macro to compute flattened array index for 3D data
#define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Compute spatial derivatives of the phase-field phi using finite differences.
 *
 * Calculates forward, central, and mixed-direction derivatives for each
 * interior grid point (excluding ghost cells) in 2D or 3D.
 *
 * @param phi     Input phase-field array of size NX*NY*NZ
 * @param fb      FieldBuffers struct holding derivative arrays to populate
 * @param params  Simulation parameters containing grid dimensions and spacings
 * @param r       Inverse grid spacings: [1/dx, 1/dy, 1/dz]
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1]
 */
void computeGradientPhi(double *phi, FieldBuffers *fb, const SimParams *params, double r[], int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    // Define interior k-range: skip ghost layers
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);

                // Forward differences
                fb->DERX_right[idx]  = (phi[idx + strides[0]] - phi[idx]) * r[0];
                fb->DERX_left[idx]   = (phi[idx] - phi[idx - strides[0]]) * r[0];
                fb->DERY_top[idx]    = (phi[idx + strides[1]] - phi[idx]) * r[1];
                fb->DERY_bottom[idx] = (phi[idx] - phi[idx - strides[1]]) * r[1];

                // Central differences
                fb->DERX_c[idx]      = 0.5 * (phi[idx + strides[0]] - phi[idx - strides[0]]) * r[0];
                fb->DERY_c[idx]      = 0.5 * (phi[idx + strides[1]] - phi[idx - strides[1]]) * r[1];

                // Mixed-direction derivatives
                fb->DERX_top[idx]    = 0.25 * ((phi[idx + strides[0]] + phi[idx + strides[0] + strides[1]] + phi[idx] + phi[idx + strides[1]])
                                         - (phi[idx - strides[0]] + phi[idx - strides[0] + strides[1]] + phi[idx] + phi[idx + strides[1]])) * r[0];
                fb->DERX_bottom[idx] = 0.25 * ((phi[idx + strides[0]] + phi[idx + strides[0] - strides[1]] + phi[idx] + phi[idx - strides[1]])
                                         - (phi[idx - strides[0]] + phi[idx - strides[0] - strides[1]] + phi[idx] + phi[idx - strides[1]])) * r[0];
                fb->DERY_right[idx]  = 0.25 * ((phi[idx + strides[1]] + phi[idx + strides[0] + strides[1]] + phi[idx] + phi[idx + strides[0]])
                                         - (phi[idx - strides[1]] + phi[idx + strides[0] - strides[1]] + phi[idx] + phi[idx + strides[0]])) * r[1];
                fb->DERY_left[idx]   = 0.25 * ((phi[idx + strides[1]] + phi[idx - strides[0] + strides[1]] + phi[idx] + phi[idx - strides[0]])
                                         - (phi[idx - strides[1]] + phi[idx - strides[0] - strides[1]] + phi[idx] + phi[idx - strides[0]])) * r[1];
            }
        }
    }
}

#undef IDX


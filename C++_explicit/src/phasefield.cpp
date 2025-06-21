#include "header.hpp"
#include <cstdlib>
#include <cmath>

// Macro to compute flattened array index for 3D data
#define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Update phase-field phi using anisotropic fluxes, free energy derivative, and noise.
 *
 * Computes directional fluxes (rj, lj, tj, bj), adds reaction term dF/dphi and noise,
 * and advances phi by one time step: phi_new = phi + dt * dphi/dt.
 *
 * @param phi     Input phase-field array.
 * @param fb      FieldBuffers containing derivatives, anisotropy, and output buffers.
 * @param params  Simulation parameters including grid dims, dt, tau, a.
 * @param r       Inverse grid spacings: [1/dx, 1/dy, 1/dz].
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1].
 */
void updatePhi(double *phi, FieldBuffers *fb, const SimParams *params, double r[], int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    double dt = params->dt;
    double a  = params->a;
    double tau= params->tau;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);

                // Compute anisotropic fluxes
                double rj = fb->ac_right[idx] * (fb->ac_right[idx] * fb->DERX_right[idx]
                                                 - fb->ac_p_right[idx] * fb->DERY_right[idx]);
                double lj = fb->ac_left[idx]  * (fb->ac_left[idx]  * fb->DERX_left[idx]
                                                 - fb->ac_p_left[idx]  * fb->DERY_left[idx]);
                double tj = fb->ac_top[idx]   * (fb->ac_top[idx]   * fb->DERY_top[idx]
                                                 + fb->ac_p_top[idx]   * fb->DERX_top[idx]);
                double bj = fb->ac_bottom[idx]* (fb->ac_bottom[idx]* fb->DERY_bottom[idx]
                                                 + fb->ac_p_bottom[idx]* fb->DERX_bottom[idx]);

                // Add noise term: a * (rand() - 0.5) scaled by phi(1-phi)
                double noise = a * ((double)std::rand() / RAND_MAX - 0.5);
                noise *= phi[idx] * (1.0 - phi[idx]);

                // Compute time derivative dphi/dt
                double dphi = (rj - lj) * r[0] + (tj - bj) * r[1];
                dphi += fb->dfdphi[idx];
                dphi += noise;
                dphi /= tau;
                fb->dphi_dt[idx] = dphi;

                // Update phi
                fb->phi_new[idx] = phi[idx] + dt * dphi;
            }
        }
    }
}

/**
 * @brief Copy updated phi values from src to dst for interior grid points.
 *
 * @param dst     Destination phi array to update.
 * @param src     Source phi_new array containing new phi values.
 * @param params  Simulation parameters for grid dims.
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1].
 */
void copyInterior(double *dst, double *src, const SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? params->Num_Z - 1 : 1;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                dst[idx] = src[idx];
            }
        }
    }
}

#undef IDX

#include "header.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>

// Macro to compute flattened array index for 3D data
#define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Update the temperature field over one time step.
 *
 * Applies diffusion via the Laplacian operator and couples to the
 * phase-field evolution (source term K * dphi/dt).
 *
 * @param temp    Input temperature array of size NX*NY*NZ.
 * @param fb      FieldBuffers containing dphi/dt and output temp_new.
 * @param params  Simulation parameters including grid dims, dt, K.
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1].
 * @param r2      Squared inverse grid spacings: [1/dx*dx, 1/dy*dy, 1/dz*dz].
 */
void updateTemp(double *temp, FieldBuffers *fb, const SimParams *params, int strides[], double r2[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    double dt = params->dt;
    double K  = params->K;

    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                // Diffusion term via Laplacian
                double lap = computeLaplacian(temp, idx, strides, r2, dim);
                // Coupling source from phase-field
                double dtemp_dt = lap + K * fb->dphi_dt[idx];
                // Time integration
                fb->temp_new[idx] = temp[idx] + dt * dtemp_dt;
            }
        }
    }
}

/**
 * @brief Compute the discrete Laplacian of arr at a given index using 5 point stencil (can be extended to 9 point stencil).
 *
 * Uses second-order central differences in each direction.
 *
 * @param arr     Input array.
 * @param index   Flattened index for which to compute Laplacian.
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1].
 * @param r2      Squared inverse grid spacings: [1/dx*dx, 1/dy*dy, 1/dz*dz].
 * @param dim     Number of spatial dimensions (2 or 3).
 * @return The Laplacian value at the specified index.
 */
double computeLaplacian(double *arr, int index, int strides[], double r2[], int dim) {
    // X-direction
    double lap = (arr[index + strides[0]] - 2.0 * arr[index] + arr[index - strides[0]]) * r2[0];
    // Y-direction
    lap += (arr[index + strides[1]] - 2.0 * arr[index] + arr[index - strides[1]]) * r2[1];
    // Z-direction if 3D
    if (dim == 3) {
        lap += (arr[index + strides[2]] - 2.0 * arr[index] + arr[index - strides[2]]) * r2[2];
    }
    return lap;
}

#undef IDX


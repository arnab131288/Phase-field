#include "header.hpp"
#include <cmath>

// Macro to compute flattened array index for 3D data
#define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Compute the derivative of free energy with respect to the phase-field φ.
 *
 * For each interior grid point (excluding ghost cells), computes:
 *   m = (alpha / PI) * atan(gamma * (T_e - temp))
 *   dfdphi = phi * (1 - phi) * (phi - 0.5 + m)
 *
 * @param phi     Input phase-field array of size NX*NY*NZ
 * @param dfdphi  Output array to store computed dF/dphi values
 * @param temp    Temperature field array
 * @param params  Simulation parameters containing grid dimensions and constants
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1]
 */
void computedfdphi(double *phi, double *dfdphi, double *temp, const SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    const double alpha = params->alpha;
    const double gamma = params->gamma;
    const double T_e   = params->T_e;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                // Coupling term: m = (alpha/PI) * atan(gamma * (T_e - temp))
                double m = (alpha / M_PI) * std::atan(gamma * (T_e - temp[idx]));
                // Derivative of free energy
                dfdphi[idx] = phi[idx] * (1.0 - phi[idx]) * (phi[idx] - 0.5 + m);
            }
        }
    }
}

#undef IDX


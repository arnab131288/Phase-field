#include "header.hpp"

// Macro to compute flattened array index for 3D data
#define IDX(i,j,k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Fill a cubic region within the array with a specified value.
 *
 * For each interior grid point, assigns vb->fillValue if the point lies
 * within the cube defined by vb->fill.c; otherwise sets it to
 * (1 - vb->fillValue).
 *
 * @param arr     Pointer to the data array to fill.
 * @param vb      Boundary and fill information for this variable.
 * @param params  Simulation parameters containing grid dimensions.
 * @param strides Array of strides for flattening 3D indices.
 */
void FillCube(double *arr, const VariableBoundary *vb, const SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    // Loop over interior grid (excluding ghost cells)
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                bool inside = (i >= vb->fill.c.x_start && i <= vb->fill.c.x_end) &&
                              (j >= vb->fill.c.y_start && j <= vb->fill.c.y_end) &&
                              (dim == 2 || (k >= vb->fill.c.z_start && k <= vb->fill.c.z_end));
                arr[idx] = inside ? vb->fillValue : (1.0 - vb->fillValue);
            }
        }
    }
}

/**
 * @brief Fill a spherical region within the array with a specified value.
 *
 * Computes the squared distance from each grid point to the sphere center;
 * if within radius, sets arr to vb->fillValue; else 1 - vb->fillValue.
 *
 * @param arr     Pointer to the data array to fill.
 * @param vb      Boundary and fill information for this variable.
 * @param params  Simulation parameters containing grid dimensions.
 * @param strides Array of strides for flattening 3D indices.
 */
void FillSphere(double *arr, const VariableBoundary *vb, const SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    // Loop over interior grid
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                double dx2 = (i - vb->fill.s.x_center) * (i - vb->fill.s.x_center);
                double dy2 = (j - vb->fill.s.y_center) * (j - vb->fill.s.y_center);
                double dz2 = (dim == 2) ? 0.0 : (k - vb->fill.s.z_center) * (k - vb->fill.s.z_center);
                double dist2 = dx2 + dy2 + dz2;
                double rad2 = vb->fill.s.radius * vb->fill.s.radius;
                arr[idx] = (dist2 <= rad2) ? vb->fillValue : (1.0 - vb->fillValue);
            }
        }
    }
}

/**
 * @brief Fill the entire interior domain with a constant value.
 *
 * Sets each interior grid point to vb->fillValue.
 *
 * @param arr     Pointer to the data array to fill.
 * @param vb      Boundary and fill information for this variable.
 * @param params  Simulation parameters containing grid dimensions.
 * @param strides Array of strides for flattening 3D indices.
 */
void FillConstant(double *arr, const VariableBoundary *vb, const SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    // Loop over interior grid
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                arr[idx] = vb->fillValue;
            }
        }
    }
}

#undef IDX

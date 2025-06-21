#include "header.hpp"

/**
 * @brief Apply boundary conditions to a 3D (or 2D) array slice.
 *
 * For each face (left/right, bottom/top, back/front), applies either
 * periodic or no-flux (Neumann) boundary conditions by copying from
 * the appropriate interior reference cell.
 *
 * @param arr     Pointer to the data array (flattened 3D)
 * @param params  Simulation parameters containing grid dimensions
 * @param strides Strides for flattening 3D indices: [NY*NZ, NZ, 1]
 * @param bc      FaceBoundary specifying conditions per direction
 */
void applyBoundaryConditions(double *arr, const SimParams *params, int strides[], FaceBoundary bc) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;

    // Macro to index flattened array
    #define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

    // X-direction boundaries
    for (int j = 1; j < NY - 1; ++j) {
        int kstart = (dim == 3) ? 1 : 0;
        int kend   = (dim == 3) ? NZ - 1 : 1;
        for (int k = kstart; k < kend; ++k) {
            int left_idx  = (dim == 2) ? (0 * strides[0] + j * strides[1]) : IDX(0, j, k);
            int right_idx = (dim == 2) ? ((NX - 1) * strides[0] + j * strides[1]) : IDX(NX - 1, j, k);
            int lref_idx  = (dim == 2) ? (1 * strides[0] + j * strides[1]) : IDX(1, j, k);
            int rref_idx  = (dim == 2) ? ((NX - 2) * strides[0] + j * strides[1]) : IDX(NX - 2, j, k);

            // Left face
            if (bc.left == BOUNDARY_PERIODIC)
                arr[left_idx] = arr[rref_idx];
            else if (bc.left == BOUNDARY_NOFLUX)
                arr[left_idx] = arr[lref_idx];

            // Right face
            if (bc.right == BOUNDARY_PERIODIC)
                arr[right_idx] = arr[lref_idx];
            else if (bc.right == BOUNDARY_NOFLUX)
                arr[right_idx] = arr[rref_idx];
        }
    }

    // Y-direction boundaries
    for (int i = 1; i < NX - 1; ++i) {
        int kstart = (dim == 3) ? 1 : 0;
        int kend   = (dim == 3) ? NZ - 1 : 1;
        for (int k = kstart; k < kend; ++k) {
            int bot_idx = (dim == 2) ? (i * strides[0] + 0 * strides[1]) : IDX(i, 0, k);
            int top_idx = (dim == 2) ? (i * strides[0] + (NY - 1) * strides[1]) : IDX(i, NY - 1, k);
            int bref_idx= (dim == 2) ? (i * strides[0] + 1 * strides[1]) : IDX(i, 1, k);
            int tref_idx= (dim == 2) ? (i * strides[0] + (NY - 2) * strides[1]) : IDX(i, NY - 2, k);

            // Bottom face
            if (bc.bottom == BOUNDARY_PERIODIC)
                arr[bot_idx] = arr[tref_idx];
            else if (bc.bottom == BOUNDARY_NOFLUX)
                arr[bot_idx] = arr[bref_idx];

            // Top face
            if (bc.top == BOUNDARY_PERIODIC)
                arr[top_idx] = arr[bref_idx];
            else if (bc.top == BOUNDARY_NOFLUX)
                arr[top_idx] = arr[tref_idx];
        }
    }

    // Z-direction boundaries (3D only)
    if (dim == 3) {
        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                int back_idx  = IDX(i, j, 0);
                int front_idx = IDX(i, j, NZ - 1);
                int bref_idx  = IDX(i, j, 1);
                int fref_idx  = IDX(i, j, NZ - 2);

                // Back face
                if (bc.back == BOUNDARY_PERIODIC)
                    arr[back_idx] = arr[fref_idx];
                else if (bc.back == BOUNDARY_NOFLUX)
                    arr[back_idx] = arr[bref_idx];

                // Front face
                if (bc.front == BOUNDARY_PERIODIC)
                    arr[front_idx] = arr[bref_idx];
                else if (bc.front == BOUNDARY_NOFLUX)
                    arr[front_idx] = arr[fref_idx];
            }
        }
    }

    #undef IDX
}

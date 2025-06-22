#include "header.h"

void applyBoundaryConditions(double *arr, SimParams *params, int strides[], FaceBoundary bc) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;

    #define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

    // X-direction
    for (int j = 1; j < NY - 1; j++) {
        int kstart = (dim == 3) ? 1 : 0;
        int kend   = (dim == 3) ? NZ - 1 : 1;
        for (int k = kstart; k < kend; k++) {
            int left_index  = (dim == 2) ? (0 * strides[0] + j * strides[1]) : IDX(0, j, k);
            int right_index = (dim == 2) ? ((NX - 1) * strides[0] + j * strides[1]) : IDX(NX - 1, j, k);
            int lref_index  = (dim == 2) ? (1 * strides[0] + j * strides[1]) : IDX(1, j, k);
            int rref_index  = (dim == 2) ? ((NX - 2) * strides[0] + j * strides[1]) : IDX(NX - 2, j, k);

            if (bc.left == BOUNDARY_PERIODIC)
                arr[left_index] = arr[rref_index];
            else if (bc.left == BOUNDARY_NOFLUX)
                arr[left_index] = arr[lref_index];

            if (bc.right == BOUNDARY_PERIODIC)
                arr[right_index] = arr[lref_index];
            else if (bc.right == BOUNDARY_NOFLUX)
                arr[right_index] = arr[rref_index];
        }
    }

    // Y-direction
    for (int i = 1; i < NX - 1; i++) {
        int kstart = (dim == 3) ? 1 : 0;
        int kend   = (dim == 3) ? NZ - 1 : 1;
        for (int k = kstart; k < kend; k++) {
            int bot_index  = (dim == 2) ? (i * strides[0] + 0 * strides[1]) : IDX(i, 0, k);
            int top_index  = (dim == 2) ? (i * strides[0] + (NY - 1) * strides[1]) : IDX(i, NY - 1, k);
            int bref_index = (dim == 2) ? (i * strides[0] + 1 * strides[1]) : IDX(i, 1, k);
            int tref_index = (dim == 2) ? (i * strides[0] + (NY - 2) * strides[1]) : IDX(i, NY - 2, k);

            if (bc.bottom == BOUNDARY_PERIODIC)
                arr[bot_index] = arr[tref_index];
            else if (bc.bottom == BOUNDARY_NOFLUX)
                arr[bot_index] = arr[bref_index];

            if (bc.top == BOUNDARY_PERIODIC)
                arr[top_index] = arr[bref_index];
            else if (bc.top == BOUNDARY_NOFLUX)
                arr[top_index] = arr[tref_index];
        }
    }

    // Z-direction (only for 3D)
    if (dim == 3) {
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                int back_index = IDX(i, j, 0);
                int front_index = IDX(i, j, NZ - 1);
                int bref_index = IDX(i, j, 1);
                int fref_index = IDX(i, j, NZ - 2);

                if (bc.back == BOUNDARY_PERIODIC)
                    arr[back_index] = arr[fref_index];
                else if (bc.back == BOUNDARY_NOFLUX)
                    arr[back_index] = arr[bref_index];

                if (bc.front == BOUNDARY_PERIODIC)
                    arr[front_index] = arr[bref_index];
                else if (bc.front == BOUNDARY_NOFLUX)
                    arr[front_index] = arr[fref_index];
            }
        }
    }

    #undef IDX
}

#include "header.h"

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])


void FillCube(double *arr, const VariableBoundary *vb, SimParams *params, int strides[]) {
    int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = kstart; k < kend; k++) {
                int index = IDX(i, j, k);
                if (i >= vb->fill.c.x_start && i <= vb->fill.c.x_end &&
                    j >= vb->fill.c.y_start && j <= vb->fill.c.y_end &&
                    (dim == 2 || (k >= vb->fill.c.x_start && k <= vb->fill.c.z_end))) {
                    arr[index] = vb->fillValue;
                } else {
                    arr[index] = 1.0 - vb->fillValue;
                }
            }
        }
    }
}

void FillSphere(double *arr, const VariableBoundary *vb, SimParams *params, int strides[]) {
    int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = kstart; k < kend; k++) {
                int index = IDX(i, j, k);
                double dx2 = (i - vb->fill.s.x_center) * (i - vb->fill.s.x_center);
                double dy2 = (j - vb->fill.s.y_center) * (j - vb->fill.s.y_center);
                double dz2 = (dim == 2) ? 0 : (k - vb->fill.s.z_center) * (k - vb->fill.s.z_center);
                if (dx2 + dy2 + dz2 <= vb->fill.s.radius * vb->fill.s.radius) {
                    arr[index] = vb->fillValue;
                } else {
                    arr[index] = 1.0 - vb->fillValue;
                }
            }
        }
    }
}

void FillConstant(double *arr, const VariableBoundary *vb, SimParams *params, int strides[]) {
    int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = kstart; k < kend; k++) {
                int index = IDX(i, j, k);
                arr[index] = vb->fillValue;
            }
        }
    }
}


#undef IDX
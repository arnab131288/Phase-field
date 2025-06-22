#include "header.h"

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])


void write_output_csv(const char *filename, double *arr, SimParams *params, int strides[]) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Warning: Could not open file %s for writing CSV output.\n", filename);
        return;
    }

    int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = kstart; k < kend; k++) {
                int index = IDX(i, j, k);
                if (dim == 2)
                    fprintf(fp, "%d,%d,%.8lf\n", i, j, arr[index]);
                else
                    fprintf(fp, "%d,%d,%d,%.8lf\n", i, j, k, arr[index]);
            }
        }
    }
    fclose(fp);
}

void write_output_vtk(const char *filename, double *arr, SimParams *params, int strides[]) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Could not open %s for writing VTK output.\n", filename);
        return;
    }

    int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Concentration output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    if (dim == 2) {
        fprintf(fp, "DIMENSIONS %d %d 1\n", NX - 2, NY - 2);
        fprintf(fp, "ORIGIN 0 0 0\n");
        fprintf(fp, "SPACING %g %g 1.0\n", params->dx, params->dy);
        fprintf(fp, "POINT_DATA %d\n", (NX - 2) * (NY - 2));
    } else {
        fprintf(fp, "DIMENSIONS %d %d %d\n", NX - 2, NY - 2, NZ - 2);
        fprintf(fp, "ORIGIN 0 0 0\n");
        fprintf(fp, "SPACING %g %g %g\n", params->dx, params->dy, params->dz);
        fprintf(fp, "POINT_DATA %d\n", (NX - 2) * (NY - 2) * (NZ - 2));
    }
    fprintf(fp, "SCALARS concentration float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;
    for (int k = kstart; k < kend; k++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int i = 1; i < NX - 1; i++) {
                int index = IDX(i, j, k);
                fprintf(fp, "%.8lf\n", arr[index]);
            }
        }
    }
    fclose(fp);
}

#undef IDX

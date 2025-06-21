#include "header.hpp"
#include <cstdio>
#include <cstdlib>

// Macro to compute flattened array index for 3D data
#define IDX(i, j, k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Write field data to a CSV file.
 *
 * In 2D: writes lines "i,j,value" for each interior point.
 * In 3D: writes lines "i,j,k,value".
 *
 * @param filename Path for CSV output.
 * @param arr      Data array of size NX*NY*NZ.
 * @param params   Simulation parameters for dimensions.
 * @param strides  Strides for flattening: [NY*NZ, NZ, 1].
 */
void write_output_csv(const char *filename, double *arr, const SimParams *params, int strides[]) {
    FILE *fp = std::fopen(filename, "w");
    if (!fp) {
        std::fprintf(stderr, "Warning: Could not open file %s for writing CSV output.\n", filename);
        return;
    }

    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = kstart; k < kend; ++k) {
                int idx = IDX(i, j, k);
                if (dim == 2) {
                    std::fprintf(fp, "%d,%d,%.8lf\n", i, j, arr[idx]);
                } else {
                    std::fprintf(fp, "%d,%d,%d,%.8lf\n", i, j, k, arr[idx]);
                }
            }
        }
    }
    std::fclose(fp);
}

/**
 * @brief Write field data to an ASCII VTK Structured Points file.
 *
 * Outputs header and then scalar values per point in row-major order.
 *
 * @param filename Path for VTK output.
 * @param arr      Data array of size NX*NY*NZ.
 * @param params   Simulation parameters for dimensions and spacing.
 * @param strides  Strides for flattening: [NY*NZ, NZ, 1].
 */
void write_output_vtk(const char *filename, double *arr, const SimParams *params, int strides[]) {
    FILE *fp = std::fopen(filename, "w");
    if (!fp) {
        std::fprintf(stderr, "Error: Could not open %s for writing VTK output.\n", filename);
        return;
    }

    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;

    // VTK header
    std::fprintf(fp, "# vtk DataFile Version 3.0\n");
    std::fprintf(fp, "Concentration output\n");
    std::fprintf(fp, "ASCII\n");
    std::fprintf(fp, "DATASET STRUCTURED_POINTS\n");

    if (dim == 2) {
        std::fprintf(fp, "DIMENSIONS %d %d 1\n", NX - 2, NY - 2);
        std::fprintf(fp, "ORIGIN 0 0 0\n");
        std::fprintf(fp, "SPACING %g %g 1.0\n", params->dx, params->dy);
        std::fprintf(fp, "POINT_DATA %d\n", (NX - 2) * (NY - 2));
    } else {
        std::fprintf(fp, "DIMENSIONS %d %d %d\n", NX - 2, NY - 2, NZ - 2);
        std::fprintf(fp, "ORIGIN 0 0 0\n");
        std::fprintf(fp, "SPACING %g %g %g\n", params->dx, params->dy, params->dz);
        std::fprintf(fp, "POINT_DATA %d\n", (NX - 2) * (NY - 2) * (NZ - 2));
    }
    std::fprintf(fp, "SCALARS Variable float 1\n");
    std::fprintf(fp, "LOOKUP_TABLE default\n");

    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;
    for (int k = kstart; k < kend; ++k) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int i = 1; i < NX - 1; ++i) {
                int idx = IDX(i, j, k);
                std::fprintf(fp, "%.8lf\n", arr[idx]);
            }
        }
    }

    std::fclose(fp);
}

#undef IDX


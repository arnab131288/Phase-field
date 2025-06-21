#include "header.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

// Macro to compute flattened array index for 3D data
#define IDX(i,j,k) ((i) * strides[0] + (j) * strides[1] + (k) * strides[2])

/**
 * @brief Read scalar field data from an ASCII VTK file into an array.
 *
 * Searches for the "LOOKUP_TABLE" marker, then reads one double per line
 * corresponding to each grid point in row-major order. Fills ghost layers
 * based on DIM to skip boundary cells.
 *
 * @param filename Path to the VTK file.
 * @param arr      Output array to populate (e.g., phi or temp).
 * @param params   Simulation parameters for grid sizing.
 * @param strides  Strides for flattening 3D indices.
 */
void read_input_vtk(const char *filename, double *arr, const SimParams *params, int strides[]) {
    FILE *fp = std::fopen(filename, "r");
    if (!fp) {
        std::fprintf(stderr, "Warning: Could not open VTK file %s for reading.\n", filename);
        std::exit(EXIT_FAILURE);
    }

    char line[256];
    // Skip header lines until LOOKUP_TABLE
    while (std::fgets(line, sizeof(line), fp)) {
        if (std::strstr(line, "LOOKUP_TABLE")) {
            break;
        }
    }

    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    // Read values for each interior point (k fixed at first slice)
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            int idx = IDX(i, j, kstart);
            if (std::fscanf(fp, "%lf", &arr[idx]) != 1) {
                std::fprintf(stderr, "Error reading data at (%d,%d).\n", i, j);
                std::fclose(fp);
                std::exit(EXIT_FAILURE);
            }
        }
    }

    std::fclose(fp);
}

/**
 * @brief Read scalar field data from a CSV file into an array.
 *
 * Expects lines formatted as "i,j,value" for each interior grid point.
 * Fills only the first slice (k = kstart) if 3D, or the single layer if 2D.
 *
 * @param filename Path to the CSV file.
 * @param arr      Output array to populate (e.g., phi or temp).
 * @param params   Simulation parameters for grid sizing.
 * @param strides  Strides for flattening 3D indices.
 */
void read_input_csv(const char *filename, double *arr, const SimParams *params, int strides[]) {
    FILE *fp = std::fopen(filename, "r");
    if (!fp) {
        std::fprintf(stderr, "Warning: Could not open CSV file %s for reading.", filename);
        std::exit(EXIT_FAILURE);
    }

    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    int kstart = (dim == 3) ? 1 : 0;

    // CSV format: i,j,value
    // Loop over interior grid points
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            int idx = IDX(i, j, kstart);
            int ii, jj;
            double value;
            if (std::fscanf(fp, "%d,%d,%lf", &ii, &jj, &value) != 3) {
                std::fprintf(stderr, "Error reading CSV data at grid (%d,%d).", i, j);
                std::fclose(fp);
                std::exit(EXIT_FAILURE);
            }
            arr[idx] = value;
        }
    }

    std::fclose(fp);
}

#undef IDX

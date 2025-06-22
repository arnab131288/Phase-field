#include "header.h"

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])

void read_input_vtk(const char *filename, double *phi, SimParams *params, int strides[]) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Warning: Could not open VTK file %s for reading.\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[256];

    // Skip lines until we find the LOOKUP_TABLE line
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "LOOKUP_TABLE")) {
            break;
        }
    }

    int kstart = (params->DIM == 3) ? 1 : 0;
    int kend   = (params->DIM == 3) ? params->Num_Z - 1 : 1;

    // Now read scalar values assuming ASCII data in row-major order
        for (int i = 1; i < params->Num_X - 1; i++) {
            for (int j = 1; j < params->Num_Y - 1; j++) {
                int index = IDX(i, j, 0);
                fscanf(fp, "%lf\n", &phi[index]); 
            }
        }
    
       fclose(fp);
}

void read_input_csv(const char *filename, double *phi, SimParams *params, int strides[]) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Warning: Could not open DAT file %s for writing.\n", filename);
        exit(EXIT_FAILURE);
    }
    
    int kstart = (params->DIM == 3) ? 1 : 0;
    int kend   = (params->DIM == 3) ? params->Num_Z - 1 : 1;

     for (int i = 1; i < params->Num_X - 1; i++) {
        for (int j = 1; j < params->Num_Y - 1; j++) {
            int index = IDX(i, j, 0);
            fscanf(fp, "%d,%d,%lf\n", &i, &j, &phi[index]);
        }
    }
     fclose(fp);   
}

#undef IDX





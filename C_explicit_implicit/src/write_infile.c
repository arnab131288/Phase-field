#include <stdio.h>
#include "header.h"

/**
 * Writes only those parameters that have been read into params,
 * in the same "KEY = VALUE" format as infile.in.
 * @param outfile Name of the output file (e.g., "outfile.in").
 * @param params Pointer to the SimParams struct populated by readInput.
 */
void writeParameters(const char *outfile, const SimParams *params) {
    FILE *fp = fopen(outfile, "w");
    if (!fp) {
        perror("Error opening output file");
        return;
    }

    // Core simulation parameters
    fprintf(fp, "DIM = %d\n", params->DIM);
    fprintf(fp, "Num_X = %d\n", params->Num_X);
    fprintf(fp, "Num_Y = %d\n", params->Num_Y);
    if (params->DIM == 3) {
        fprintf(fp, "Num_Z = %d\n", params->Num_Z);
    }
    fprintf(fp, "dx = %g\n", params->dx);
    fprintf(fp, "dy = %g\n", params->dy);
    if (params->DIM == 3) {
        fprintf(fp, "dz = %g\n", params->dz);
    }
    fprintf(fp, "dt = %g\n", params->dt);
    fprintf(fp, "total_steps = %d\n", params->total_timesteps);
    fprintf(fp, "timebreak = %d\n", params->timebreak);
    fprintf(fp, "epsilon = %g\n", params->epsilon);
    fprintf(fp, "tau = %g\n", params->tau);
    fprintf(fp, "delta = %g\n", params->delta);
    fprintf(fp, "j = %d\n", params->j);
    fprintf(fp, "alpha = %g\n", params->alpha);
    fprintf(fp, "gamma = %g\n", params->gamma);
    fprintf(fp, "a = %g\n", params->a);
    fprintf(fp, "K = %g\n", params->K);
    fprintf(fp, "T_e = %g\n", params->T_e);

    for (int i = 0; i < params->numVariables; ++i) {
         const VariableBoundary *vb = &params->variables[i];
         if (vb->fillType == FILL_CUBE) {
             // Format: Fill_Cube = varName,x0,y0,x1,y1[,z0,z1];
             if (params->DIM == 2) {
                 fprintf(fp,
                     "Fill_Cube = %s,%d,%d,%d,%d;\n",
                     vb->varName,
                     vb->fill.c.x_start,
                     vb->fill.c.y_start,
                     vb->fill.c.x_end,
                     vb->fill.c.y_end
                 );
             } else {
                 fprintf(fp,
                     "Fill_Cube = %s,%d,%d,%d,%d,%d,%d;\n",
                     vb->varName,
                     vb->fill.c.x_start,
                     vb->fill.c.y_start,
                     vb->fill.c.x_end,
                     vb->fill.c.y_end,
                     vb->fill.c.z_start,
                     vb->fill.c.z_end
                 );
             }
         }
         else if (vb->fillType == FILL_SPHERE) {
             // Format: Fill_Sphere = varName,radius,x_center,y_center[,z_center];
             if (params->DIM == 2) {
                 fprintf(fp,
                     "Fill_Sphere = %s,%d,%d,%d;\n",
                     vb->varName,
                     vb->fill.s.radius,
                     vb->fill.s.x_center,
                     vb->fill.s.y_center
                 );
             } else {
                 fprintf(fp,
                     "Fill_Sphere = %s,%d,%d,%d,%d;\n",
                     vb->varName,
                     vb->fill.s.radius,
                     vb->fill.s.x_center,
                     vb->fill.s.y_center,
                     vb->fill.s.z_center
                 );
             }
         }
    }

    // Boundary condition mode
    // Boundary conditions for each variable
    const char *bcNames[] = {"UNDEFINED","NOFLUX","PERIODIC"};
    for (int i = 0; i < params->numVariables; i++) {
        const VariableBoundary *vb = &params->variables[i];
        FaceBoundary bc = vb->bc;
        // e.g. “boundary = conc,NOFLUX,PERIODIC,NOFLUX,NOFLUX” (+ front/back if 3D)
        fprintf(fp, "boundary = %s,%s,%s,%s,%s",
            vb->varName,
            bcNames[bc.top],
            bcNames[bc.bottom],
            bcNames[bc.left],
            bcNames[bc.right]);
        if (params->DIM == 3) {
            fprintf(fp, ",%s,%s", bcNames[bc.front], bcNames[bc.back]);
        }
        fprintf(fp, "\n");
    }

    if (params->RESPAWN) fprintf(fp, "RESPAWN = %d\n", params->RESPAWN);
    if (params->RESPAWN) {
        fprintf(fp, "restart_time = %d\n", params->restart_time);
    }

    if (params->WRITE_TO_VTK) fprintf(fp, "WRITE_TO_VTK = %d\n", params->WRITE_TO_VTK);
    if (params->WRITE_TO_CSV) fprintf(fp, "WRITE_TO_CSV = %d\n", params->WRITE_TO_CSV);

    // Add any other SimParams fields here similarly...

    fclose(fp);
}

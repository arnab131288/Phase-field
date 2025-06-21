#include "header.hpp"
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>

/**
 * @brief Write simulation parameters back to an output file in "KEY = VALUE" format.
 *
 * Only parameters that were originally read into params are written.
 * Output follows the same syntax as the input file (infile.in).
 *
 * @param outfile Path to the output file (e.g., "outfile.in").
 * @param params  Pointer to populated SimParams struct.
 */
void writeParameters(const char *outfile, const SimParams *params) {
    // Open file for writing
    FILE *fp = std::fopen(outfile, "w");
    if (!fp) {
        std::fprintf(stderr, "Error opening output file '%s': %s\n", outfile, std::strerror(errno));
        return;
    }

    // Core simulation parameters
    std::fprintf(fp, "DIM = %d\n", params->DIM);
    std::fprintf(fp, "Num_X = %d\n", params->Num_X);
    std::fprintf(fp, "Num_Y = %d\n", params->Num_Y);
    if (params->DIM == 3) {
        std::fprintf(fp, "Num_Z = %d\n", params->Num_Z);
    }
    std::fprintf(fp, "dx = %g\n", params->dx);
    std::fprintf(fp, "dy = %g\n", params->dy);
    if (params->DIM == 3) {
        std::fprintf(fp, "dz = %g\n", params->dz);
    }
    std::fprintf(fp, "dt = %g\n", params->dt);
    std::fprintf(fp, "total_steps = %d\n", params->total_timesteps);
    std::fprintf(fp, "timebreak = %d\n", params->timebreak);
    std::fprintf(fp, "epsilon = %g\n", params->epsilon);
    std::fprintf(fp, "tau = %g\n", params->tau);
    std::fprintf(fp, "delta = %g\n", params->delta);
    std::fprintf(fp, "j = %d\n", params->j);
    std::fprintf(fp, "alpha = %g\n", params->alpha);
    std::fprintf(fp, "gamma = %g\n", params->gamma);
    std::fprintf(fp, "a = %g\n", params->a);
    std::fprintf(fp, "K = %g\n", params->K);
    std::fprintf(fp, "T_e = %g\n", params->T_e);

    // Fill definitions
    for (int i = 0; i < params->numVariables; ++i) {
        const VariableBoundary &vb = params->variables[i];
        if (vb.fillType == FILL_CUBE) {
            // Fill_Cube = varName,fval,x0,x1,y0,y1[,z0,z1]
            if (params->DIM == 2) {
                std::fprintf(fp,
                             "Fill_Cube = %s,%g,%d,%d,%d,%d;\n",
                             vb.varName,
                             vb.fillValue,
                             vb.fill.c.x_start,
                             vb.fill.c.x_end,
                             vb.fill.c.y_start,
                             vb.fill.c.y_end);
            } else {
                std::fprintf(fp,
                             "Fill_Cube = %s,%g,%d,%d,%d,%d,%d,%d;\n",
                             vb.varName,
                             vb.fillValue,
                             vb.fill.c.x_start,
                             vb.fill.c.x_end,
                             vb.fill.c.y_start,
                             vb.fill.c.y_end,
                             vb.fill.c.z_start,
                             vb.fill.c.z_end);
            }
        } else if (vb.fillType == FILL_SPHERE) {
            // Fill_Sphere = varName,fval,radius,xc,yc[,zc]
            if (params->DIM == 2) {
                std::fprintf(fp,
                             "Fill_Sphere = %s,%g,%d,%d,%d;\n",
                             vb.varName,
                             vb.fillValue,
                             vb.fill.s.radius,
                             vb.fill.s.x_center,
                             vb.fill.s.y_center);
            } else {
                std::fprintf(fp,
                             "Fill_Sphere = %s,%g,%d,%d,%d,%d;\n",
                             vb.varName,
                             vb.fillValue,
                             vb.fill.s.radius,
                             vb.fill.s.x_center,
                             vb.fill.s.y_center,
                             vb.fill.s.z_center);
            }
        } else if (vb.fillType == FILL_CONSTANT) {
            // Fill_Constant = varName,fval
            std::fprintf(fp,
                         "Fill_Constant = %s,%g;\n",
                         vb.varName,
                         vb.fillValue);
        }
    }

    // Boundary conditions
    const char *bcNames[] = {"UNDEFINED", "NOFLUX", "PERIODIC"};
    for (int i = 0; i < params->numVariables; ++i) {
        const VariableBoundary &vb = params->variables[i];
        const FaceBoundary &bc = vb.bc;
        std::fprintf(fp,
                     "boundary = %s,%s,%s,%s,%s",
                     vb.varName,
                     bcNames[bc.top],
                     bcNames[bc.bottom],
                     bcNames[bc.left],
                     bcNames[bc.right]);
        if (params->DIM == 3) {
            std::fprintf(fp, ",%s,%s",
                         bcNames[bc.front],
                         bcNames[bc.back]);
        }
        std::fprintf(fp, "\n");
    }

    // Respawn and output options
    if (params->RESPAWN) {
        std::fprintf(fp, "RESPAWN = %d\n", params->RESPAWN);
        std::fprintf(fp, "restart_time = %d\n", params->restart_time);
    }
    if (params->WRITE_TO_VTK)   std::fprintf(fp, "WRITE_TO_VTK = %d\n", params->WRITE_TO_VTK);
    if (params->WRITE_TO_CSV)   std::fprintf(fp, "WRITE_TO_CSV = %d\n", params->WRITE_TO_CSV);

    // Close file
    std::fclose(fp);
}

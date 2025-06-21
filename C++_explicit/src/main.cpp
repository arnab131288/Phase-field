#include "header.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <unistd.h>

/*
 * main.cpp
 *
 * Entry point for the simulation application. Performs the following:
 *  - Parses command-line arguments for an input configuration file
 *  - Reads simulation parameters from the input file
 *  - Manages output directory creation and optional cleanup
 *  - Initializes simulation variables and field buffers
 *  - Handles respawn logic: loading previous phi and temperature fields
 *  - Executes the main time-stepping loop:
 *      a) Applies boundary conditions
 *      b) Computes free energy derivatives
 *      c) Computes gradients and anisotropy
 *      d) Updates phase-field and temperature fields
 *      e) Periodically writes output in VTK or CSV formats
 *  - Cleans up allocated memory on exit
 */

int main(int argc, char* argv[]) {
    // Validate command-line arguments
    if (argc < 2) {
        std::fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Initialize simulation parameters
    SimParams params{};
    params.numVariables = 0;
    if (readParameters(argv[1], &params) != 0) {
        // Exit if parameter reading fails
        return EXIT_FAILURE;
    }

    const char* folder = "output";

    // Handle output directory creation or cleanup
    if (!params.RESPAWN) {
        // Check if folder exists
        if (access(folder, F_OK) != -1) {
            char response[4];
            std::printf("Folder '%s' already exists. Delete it? (yes/no): ", folder);
            int matches = std::scanf("%3s", response);
            if (matches != 1) {
                std::fprintf(stderr, "Input failed or EOF\n");
                // handle error or exit:
                return 1;
            }

            if (std::strcmp(response, "yes") == 0) {
                // Remove existing directory
                char command[256];
                std::snprintf(command, sizeof(command), "rm -rf %s", folder);
                if (std::system(command) != 0) {
                    std::perror("Failed to delete folder");
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::printf("Folder not deleted. Exiting.\n");
                std::exit(EXIT_SUCCESS);
            }
        }
        // Create the output folder
        if (mkdir(folder, 0777) == -1) {
            std::perror("Error creating folder");
            std::exit(EXIT_FAILURE);
        }
        std::printf("Folder '%s' created successfully.\n", folder);
    }

    // Save updated parameters
    writeParameters("output/outfile.in", &params);

    // Precompute inverse grid spacing and squared spacing
    double r[MAX_DIM] = { 1.0 / params.dx, 1.0 / params.dy, 1.0 / params.dz };
    double r2[MAX_DIM] = { 1.0/(params.dx*params.dx), 1.0/(params.dy*params.dy), 1.0/(params.dz*params.dz) };
    int strides[MAX_DIM] = { params.Num_Y * params.Num_Z, params.Num_Z, 1 };

    // Allocate and register variable arrays
    setupVariables(&params);

    // Retrieve primary field arrays
    double* phi = getDataArray("phi");
    if (!phi) {
        std::fprintf(stderr, "Error: No data array for 'phi'.\n");
        return EXIT_FAILURE;
    }
    double* temp = getDataArray("temp");
    if (!temp) {
        std::fprintf(stderr, "Error: No data array for 'temp'.\n");
        return EXIT_FAILURE;
    }

    // Allocate buffers for intermediate computations
    FieldBuffers fb;
    allocateFieldBuffers(&params, &fb);

    // Initial condition: fill or respawn fields
    if (!params.RESPAWN) {
        // Fill fields based on defined shapes
        for (int i = 0; i < params.numVariables; ++i) {
            VariableBoundary* vb = &params.variables[i];
            double* arr = getDataArray(vb->varName);
            switch (vb->fillType) {
                case FILL_CUBE:    FillCube(arr, vb, &params, strides); break;
                case FILL_SPHERE:  FillSphere(arr, vb, &params, strides); break;
                case FILL_CONSTANT:FillConstant(arr, vb, &params, strides); break;
                default: break;
            }
        }
    } else {
        // Respawn: read previous fields from VTK or CSV
        char filename[256];
        if (params.WRITE_TO_VTK) {
            std::fprintf(stderr, "Reading phi from output/phi.%d.vtk\n", params.restart_time);
            std::snprintf(filename, sizeof(filename), "output/phi_%d.vtk", params.restart_time);
            read_input_vtk(filename, phi, &params, strides);
            std::fprintf(stderr, "Reading temp from output/temp.%d.vtk\n", params.restart_time);
            std::snprintf(filename, sizeof(filename), "output/temp_%d.vtk", params.restart_time);
            read_input_vtk(filename, temp, &params, strides);
        } else if (params.WRITE_TO_CSV) {
            std::fprintf(stderr, "Reading phi from output/phi.%d.csv\n", params.restart_time);
            std::snprintf(filename, sizeof(filename), "output/phi_%d.csv", params.restart_time);
            read_input_csv(filename, phi, &params, strides);
            std::fprintf(stderr, "Reading temp from output/temp.%d.csv\n", params.restart_time);
            std::snprintf(filename, sizeof(filename), "output/temp_%d.csv", params.restart_time);
            read_input_csv(filename, temp, &params, strides);
        }
    }

    // Write initial output if not respawning
    if (!params.RESPAWN) {
        if (params.WRITE_TO_VTK) {
            write_output_vtk("output/phi_0.vtk", phi, &params, strides);
            write_output_vtk("output/temp_0.vtk", temp, &params, strides);
        } else if (params.WRITE_TO_CSV) {
            write_output_csv("output/phi_0.csv", phi, &params, strides);
            write_output_csv("output/temp_0.csv", temp, &params, strides);
        }
    }

    // Main simulation loop over timesteps
    for (int t = 1; t <= params.total_timesteps; ++t) {
        // a) Apply boundary conditions to phi
        if (auto vb_phi = findVariableBoundary("phi", &params)) {
            applyBoundaryConditions(getDataArray("phi"), &params, strides, vb_phi->bc);
        }
        // b) Compute free-energy derivative
        computedfdphi(phi, fb.dfdphi, temp, &params, strides);
        // c) Compute gradients and anisotropy
        computeGradientPhi(phi, &fb, &params, r, strides);
        computeAnisotropy(&fb, &params, strides);
        // d) Update phi
        updatePhi(phi, &fb, &params, r, strides);
        // e) Apply boundary conditions to temp
        if (auto vb_temp = findVariableBoundary("temp", &params)) {
            applyBoundaryConditions(getDataArray("temp"), &params, strides, vb_temp->bc);
        }
        // f) Update temp
        updateTemp(temp, &fb, &params, strides, r2);
        // g) Copy new values back to main arrays
        copyInterior(phi, fb.phi_new, &params, strides);
        copyInterior(temp, fb.temp_new, &params, strides);
        // h) Periodic output
        if (t % params.timebreak == 0) {
            int t0 = params.RESPAWN ? params.restart_time : 0;
            char filename[256];
            if (params.WRITE_TO_VTK) {
                std::snprintf(filename, sizeof(filename), "output/phi_%d.vtk", t + t0);
                write_output_vtk(filename, phi, &params, strides);
                std::snprintf(filename, sizeof(filename), "output/temp_%d.vtk", t + t0);
                write_output_vtk(filename, temp, &params, strides);
                std::printf("Step %d: VTK output complete\n", t + t0);
            } else if (params.WRITE_TO_CSV) {
                std::snprintf(filename, sizeof(filename), "output/phi_%d.csv", t + t0);
                write_output_csv(filename, phi, &params, strides);
                std::snprintf(filename, sizeof(filename), "output/temp_%d.csv", t + t0);
                write_output_csv(filename, temp, &params, strides);
                std::printf("Step %d: CSV output complete\n", t + t0);
            }
        }
    }

    // Cleanup allocated memory and exit
    freeGlobalVariableArrays();
    freeFieldBuffers(&fb);
    return EXIT_SUCCESS;
}

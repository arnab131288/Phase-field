#include "header.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    
    SimParams params = {0};
    params.numVariables = 0;
    if (readParameters(argv[1], &params) != 0) {
        return EXIT_FAILURE;
    }


    const char *folder = "output";
    
    if(!params.RESPAWN) {
        // Check if folder exists
        if (access(folder, F_OK) != -1) {
        // Folder exists
        char response[4];
        printf("Folder '%s' already exists. Do you want to delete it? (yes/no): ", folder);
        scanf("%3s", response);

        if (strcmp(response, "yes") == 0) {
            // Remove the directory and its contents
            char command[256];
            snprintf(command, sizeof(command), "rm -rf %s", folder);
            int ret = system(command);
            if (ret != 0) {
                perror("Failed to delete folder");
                exit(EXIT_FAILURE);
            }
        } else {
            printf("Folder not deleted. Exiting.\n");
            exit(0);
          }
       }

         // Create the folder
         if (mkdir(folder, 0777) == -1) {
             perror("Error creating folder");
             exit(0);
        }

        printf("Folder '%s' created successfully.\n", folder);
    }
    

    writeParameters("output/outfile.in", &params);

    double r[MAX_DIM] = {
        1.0 / (params.dx),
        1.0 / (params.dy),
        1.0 / (params.dz)
    };
    
    double r2[MAX_DIM] = {
        1.0 / (params.dx * params.dx),
        1.0 / (params.dy * params.dy),
        1.0 / (params.dz * params.dz)
    };

    int strides[MAX_DIM] = {
        params.Num_Y * params.Num_Z,
        params.Num_Z,
        1
    };
    

    
    // Automatically allocate and register each variable based on boundary definitions.
    setupVariables(&params);

    double *phi = getDataArray("phi");
    if (!phi) {
        fprintf(stderr, "Error: No data array for variable 'phi'.\n");
        return EXIT_FAILURE;
    }

    double *temp = getDataArray("temp");
    if (!phi) {
        fprintf(stderr, "Error: No data array for variable 'temp'.\n");
        return EXIT_FAILURE;
    }

    FieldBuffers fb;
    allocateFieldBuffers(&params, &fb);


    
    if (!params.RESPAWN) {
        for (int i = 0; i < params.numVariables; i++) {
            VariableBoundary *vb = &params.variables[i];
            if (vb->fillType == FILL_CUBE) {
               double *arr = getDataArray(vb->varName);
               FillCube(arr, vb, &params, strides);
            } else if (vb->fillType == FILL_SPHERE) {
                double *arr = getDataArray(vb->varName);
                FillSphere(arr, vb, &params, strides);
            } else if (vb->fillType == FILL_CONSTANT) {
                double *arr = getDataArray(vb->varName);
                FillConstant(arr, vb, &params, strides);
            }
       }
     } else {
           char filename[256];
           if (params.WRITE_TO_VTK) {
              fprintf(stderr, "Reading input from VTK file output/phi.%d.vtk\n",params.restart_time);
              snprintf(filename, sizeof(filename), "output/phi_%d.vtk", params.restart_time);
              read_input_vtk(filename, phi, &params, strides);
              fprintf(stderr, "Reading input from VTK file output/temp.%d.vtk\n",params.restart_time);
              snprintf(filename, sizeof(filename), "output/temp_%d.vtk", params.restart_time);
              read_input_vtk(filename, temp, &params, strides);
           } else if (params.WRITE_TO_CSV) {
               fprintf(stderr, "Reading input from CSV file output/phi.%d.csv\n",params.restart_time);
               snprintf(filename, sizeof(filename), "output/phi_%d.csv", params.restart_time);
               read_input_csv(filename, phi, &params, strides);
               fprintf(stderr, "Reading input from CSV file output/temp.%d.csv\n",params.restart_time);
               snprintf(filename, sizeof(filename), "output/temp_%d.csv", params.restart_time);
               read_input_csv(filename, temp, &params, strides);
           }         
      }
    
    

    if (!params.RESPAWN) {
        if (params.WRITE_TO_VTK) {
            write_output_vtk("output/phi_0.vtk", phi, &params, strides);
            write_output_vtk("output/temp_0.vtk", temp, &params, strides);
        } else if (params.WRITE_TO_CSV) {
            write_output_csv("output/phi_0.csv", phi, &params, strides);
            write_output_csv("output/temp_0.csv", temp, &params, strides);
        }
    }
   

    /* Main time-stepping loop */
    for (int t = 1; t <= params.total_timesteps; t++) {
        
        // Step a) Apply BCs to phase-field
        VariableBoundary *vb_phi = findVariableBoundary("phi", &params);
        if (vb_phi) {
           double *dataArray = getDataArray("phi");
           applyBoundaryConditions(dataArray, &params, strides, vb_phi->bc);
        }
        
        /* b) Compute derivative of free energy */
        computedfdphi(phi, fb.dfdphi, temp, &params, strides);
        
        /* c) Compute gradients of phase-field and anisotropy function*/
        computeGradientPhi(phi, &fb, &params, r, strides);
        
        computeAnisotropy(&fb, &params, strides);
         
        /* d) Update phase-field */
        updatePhi(phi, &fb, &params, r, strides);

        /* e) Apply BCs to temperature field*/
         VariableBoundary *vb_temp = findVariableBoundary("temp", &params);
         if (vb_temp) {
             double *dataArray = getDataArray("temp");
           applyBoundaryConditions(dataArray, &params, strides, vb_temp->bc);
         }

        /* f) Update temperature field*/
        updateTemp(temp, &fb, &params, strides, r2); 


        /* g) Copy phi_new back to phi */
        copyInterior(phi, fb.phi_new, &params, strides);
        /* h) Copy temp_new back to temp */
        copyInterior(temp, fb.temp_new, &params, strides);

        
        /* f) Periodic output */
        if (t % params.timebreak == 0) {
            char filename[256];
            int t0 = 0;
            if (params.RESPAWN) {
                t0 = params.restart_time;
            }
             if (params.WRITE_TO_VTK) {
                snprintf(filename, sizeof(filename), "output/phi_%d.vtk", t + t0);
                write_output_vtk(filename, phi, &params, strides);
                snprintf(filename, sizeof(filename), "output/temp_%d.vtk", t + t0);
                write_output_vtk(filename, temp, &params, strides);
                printf("Step %d: file writing complete\n", t + t0);

            } else if (params.WRITE_TO_CSV) {
                snprintf(filename, sizeof(filename), "output/phi_%d.csv", t + t0);
                write_output_csv(filename, phi, &params, strides);
                snprintf(filename, sizeof(filename), "output/temp_%d.csv", t + t0);
                write_output_csv(filename, temp, &params, strides);
                printf("Step %d: file writing complete\n", t + t0);
            }
            
        }
    }

    /* Free memory */
    freeGlobalVariableArrays();
    freeFieldBuffers(&fb);

    return EXIT_SUCCESS;
}



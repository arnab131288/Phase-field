#include "header.h"


// Helper function to trim leading and trailing whitespace in a string.
void trim(char *str) {
    char *start = str;
    char *end;

    // Trim leading whitespace.
    while (*start && isspace((unsigned char)*start)) {
        start++;
    }

    if (start != str) {
        memmove(str, start, strlen(start) + 1);
    }
    
    // Trim trailing whitespace.
    end = str + strlen(str) - 1;
    while (end >= str && isspace((unsigned char)*end)) {
        *end = '\0';
        end--;
    }
}

int readParameters(const char *filename, SimParams *params) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file '%s'\n", filename);
        return 1;
    }

    /* Flags to check that each parameter is read */
    int found_DIM = 0, found_Num_X = 0, found_Num_Y = 0, found_Num_Z = 0, found_dx = 0, found_dy = 0, found_dz = 0;
    int found_dt = 0, found_total_steps = 0, found_timebreak = 0;
    int found_epsilon = 0, found_tau = 0, found_delta = 0, found_j = 0, found_theta_0 = 0;
    int found_alpha = 0, found_gamma = 0, found_a = 0, found_K = 0, found_T_e = 0, found_boundary = 0;
    int found_fill_cube = 0, found_fill_sphere = 0, found_fill_constant = 0;
    int found_RESPAWN = 0, found_restarttime = 0;
    int found_write_to_csv = 0, found_write_to_dat = 0, found_write_to_vtk = 0;

    char line[256];
    
    while (fgets(line, sizeof(line), fp)) {
        /* Skip comment lines (lines starting with '#' or empty lines) */
        if (line[0] == '#' || line[0] == '\n')
            continue;
        
        // Parse key and value separated by '=' and ending with a ';'
        char key[64], value[64];
        if (!strchr(line, ';') || sscanf(line, " %63[^=]=%63[^;];", key, value) != 2) {
            fprintf(stderr, "Warning: Could not parse line (missing '=' or ';'?): %s", line);
            exit(EXIT_FAILURE);
        }


        
        // Trim whitespace from key and value.
        trim(key);
        trim(value);

        // For debugging, uncomment the next line.
        // printf("DEBUG: key='%s', value='%s'\n", key, value);

        // For each possible key, convert the value string accordingly.
        if (strcmp(key, "DIM") == 0) {
            params->DIM = atoi(value);
            found_DIM = 1;
        } else if (strcmp(key, "Num_X") == 0) {
            params->Num_X = atoi(value);
            found_Num_X = 1;
        } else if (strcmp(key, "Num_Y") == 0) {
            params->Num_Y = atoi(value);
            found_Num_Y = 1;
        } else if (strcmp(key, "Num_Z") == 0) {
            params->Num_Z = atoi(value);
            found_Num_Z = 1;
        } else if (strcmp(key, "dx") == 0) {
            params->dx = atof(value);
            found_dx = 1;
        } else if (strcmp(key, "dy") == 0) {
            params->dy = atof(value);
            found_dy = 1;
        } else if (strcmp(key, "dz") == 0) {
            params->dz = atof(value);
            found_dz = 1;
        } else if (strcmp(key, "dt") == 0) {
            params->dt = atof(value);
            found_dt = 1;
        } else if (strcmp(key, "total_steps") == 0) {
            params->total_timesteps = atoi(value);
            found_total_steps = 1;
        } else if (strcmp(key, "timebreak") == 0) {
            params->timebreak = atoi(value);
            found_timebreak = 1;
        } else if (strcmp(key, "epsilon") == 0) {
            params->epsilon = atof(value);
            found_epsilon = 1;
        } else if (strcmp(key, "tau") == 0) {
            params->tau = atof(value);
            found_tau = 1;
        } else if (strcmp(key, "delta") == 0) {
            params->delta = atof(value);
            found_delta = 1;
        } else if (strcmp(key, "j") == 0) {
            params->j = atof(value);
            found_j = 1;
        } else if (strcmp(key, "theta_0") == 0) {
            params->theta_0 = atof(value);
            found_theta_0 = 1;
        } else if (strcmp(key, "alpha") == 0) {
            params->alpha = atof(value);
            found_alpha = 1;
        } else if (strcmp(key, "gamma") == 0) {
            params->gamma = atof(value);
            found_gamma = 1;
        } else if (strcmp(key, "a") == 0) {
            params->a = atof(value);
            found_a = 1;
        } else if (strcmp(key, "K") == 0) {
            params->K = atof(value);
            found_K = 1;
        } else if (strcmp(key, "T_e") == 0) {
            params->T_e = atof(value);
            found_T_e = 1;
        } else if (strcmp(key, "boundary") == 0) {
            // Expected format: variable,top,bottom,left,right
            char *token;
            int tokenCount = 0;
            char varName[MAX_VAR_NAME] = "";
            
            // Temporary FaceBoundary to hold the parsed boundary conditions.
            FaceBoundary tempBC;
            // Initialize to undefined (optional, but helps with debugging).
            tempBC.top = tempBC.bottom = tempBC.left = tempBC.right = tempBC.front = tempBC.back = BOUNDARY_UNDEFINED;
            
            token = strtok(value, ",");
            while (token != NULL) {
                trim(token);
                tokenCount++;
                if (tokenCount == 1) {
                    // First token is the variable name.
                    strncpy(varName, token, MAX_VAR_NAME - 1);
                    varName[MAX_VAR_NAME - 1] = '\0';
                } else {
                    BoundaryType bc;
                    if (strcmp(token, "NOFLUX") == 0)
                        bc = BOUNDARY_NOFLUX;
                    else if (strcmp(token, "PERIODIC") == 0)
                        bc = BOUNDARY_PERIODIC;
                    else {
                        fprintf(stderr, "Error: unknown boundary condition '%s' for token %d.\n", token, tokenCount);
                        fclose(fp);
                        return 1;
                    }
                    
                   switch (tokenCount) {
                        case 2: tempBC.top = bc; break;
                        case 3: tempBC.bottom = bc; break;
                        case 4: tempBC.left = bc; break;
                        case 5: tempBC.right = bc; break;
                        case 6: if (params->DIM == 3) tempBC.front = bc; break;
                        case 7: if (params->DIM == 3) tempBC.back = bc; break;
                        default:
                            fprintf(stderr, "Warning: Extra token '%s' in boundary specification for variable '%s'; ignoring.\n",
                                    token, varName);
                    }
                }
                token = strtok(NULL, ",");
            }

            if ((params->DIM == 2 && tokenCount != 5) || (params->DIM == 3 && tokenCount != 7)) {
                fprintf(stderr, "Error: Expected %d tokens for 'boundary' but got %d.\n",
                        (params->DIM == 2) ? 5 : 7, tokenCount);
                fclose(fp);
                return 1;
            }
            
            
            // Now update the SimParams structure.
            // Check if an entry for this variable already exists.
            int found = 0;
            for (int i = 0; i < params->numVariables; i++) {
                if (strcmp(params->variables[i].varName, varName) == 0) {
                    // Update boundary conditions for this variable.
                    params->variables[i].bc = tempBC;
                    found = 1;
                    break;
                }
            }
            // If no entry exists, add a new one.
            if (!found) {
                if (params->numVariables >= MAX_VARIABLES) {
                    fprintf(stderr, "Error: Exceeded maximum number of variables (%d).\n", MAX_VARIABLES);
                    fclose(fp);
                    return 1;
                }
                strncpy(params->variables[params->numVariables].varName, varName, MAX_VAR_NAME - 1);
                params->variables[params->numVariables].varName[MAX_VAR_NAME - 1] = '\0';
                params->variables[params->numVariables].bc = tempBC;
                params->numVariables++;
            }
            found_boundary = 1;
        

        /*Fill-Cube parser (expects 6 tokens in 2D or 8 in 3D)*/
   
        } else if (strcasecmp(key, "Fill_Cube") == 0) {
            // 1) count tokens
            int commaCount = 0;
            for (char *p = value; *p; ++p) if (*p == ',') ++commaCount;
                 int tokenCount = commaCount + 1;
                 int expected = (params->DIM == 3 ? 8 : 6);
                 if (tokenCount != expected) {
                    fprintf(stderr, "Error: Fill_Cube expected %d tokens but got %d.\n", expected, tokenCount);
                    fclose(fp);
                    return 1;
                 }

            // 2) parse them
            char  *tok, *save = NULL;
            int    n = 0;
            char   name[MAX_VAR_NAME];
            double fval = 0.0;
            int    x0=0, x1=0, y0=0, y1=0, z0=0, z1=0;
            for (tok = strtok_r(value, ",", &save);
                 tok;
                 tok = strtok_r(NULL, ",", &save), ++n)
            {
                 trim(tok);
                 switch (n) {
                      case 0: strncpy(name, tok, MAX_VAR_NAME-1); name[MAX_VAR_NAME-1] = '\0'; break;
                      case 1: fval = atof(tok); break;
                      case 2: x0   = atoi(tok); break;
                      case 3: x1   = atoi(tok); break;
                      case 4: y0   = atoi(tok); break;
                      case 5: y1   = atoi(tok); break;
                      case 6: if (params->DIM==3) z0 = atoi(tok); break;
                      case 7: if (params->DIM==3) z1 = atoi(tok); break;
                      default: break;
                    }
            }

            // 3) lookup & assign
            VariableBoundary *vb = findVariableBoundary(name, params);
            if (!vb) {
                fprintf(stderr, "Error: Fill_Cube for undeclared '%s'.\n", name);
                fclose(fp);
                return 1;
            }
            vb->fillType       = FILL_CUBE;
            vb->fillValue      = fval;
            vb->fill.c.x_start = x0;
            vb->fill.c.x_end   = x1;
            vb->fill.c.y_start = y0;
            vb->fill.c.y_end   = y1;
            vb->fill.c.z_start = z0;
            vb->fill.c.z_end   = z1;
            found_fill_cube = 1;
        


   /*Fill-Sphere parser (expects 5 tokens in 2D; 6 in 3D)*/
   
        } else if (strcasecmp(key, "Fill_Sphere") == 0) {
            // 1) count tokens
            int commaCount = 0;
            for (char *p = value; *p; ++p) if (*p == ',') ++commaCount;
            int tokenCount = commaCount + 1;
            int expected = (params->DIM == 3 ? 6 : 5);
            if (tokenCount != expected) {
                fprintf(stderr, "Error: Fill_Sphere expected %d tokens but got %d.\n", expected, tokenCount);
                fclose(fp);
                return 1;
            }

            // 2) parse
            char  *tok, *save = NULL;
            int    n = 0;
            char   name[MAX_VAR_NAME];
            double fval = 0.0, radius = 0.0;
            int    cx=0, cy=0, cz=0;
            for (tok = strtok_r(value, ",", &save);
                 tok;
                 tok = strtok_r(NULL, ",", &save), ++n)
            {
                 trim(tok);
                 switch (n) {
                      case 0: strncpy(name, tok, MAX_VAR_NAME-1); name[MAX_VAR_NAME-1] = '\0'; break;
                      case 1: fval   = atof(tok); break;
                      case 2: radius = atof(tok); break;
                      case 3: cx     = atoi(tok); break;
                      case 4: cy     = atoi(tok); break;
                      case 5: if (params->DIM==3) cz = atoi(tok); break;
                      default: break;
                    }
            }

            // 3) lookup & assign
            VariableBoundary *vb = findVariableBoundary(name, params);
            if (!vb) {
                fprintf(stderr, "Error: Fill_Sphere for undeclared '%s'.\n", name);
                fclose(fp);
                return 1;
            }
            vb->fillType        = FILL_SPHERE;
            vb->fillValue       = fval;
            vb->fill.s.x_center = cx;
            vb->fill.s.y_center = cy;
            vb->fill.s.z_center = cz;
            vb->fill.s.radius   = (int)radius;
            found_fill_sphere   = 1;




   /*Fill-Constant parser (always expects exactly 2 tokens)*/
  
        } else if (strcasecmp(key, "Fill_Constant") == 0) {
            // 1) count tokens
            int commaCount = 0;
            for (char *p = value; *p; ++p) if (*p == ',') ++commaCount;
            int tokenCount = commaCount + 1;
            if (tokenCount != 2) {
                fprintf(stderr, "Error: Fill_Constant expected 2 tokens but got %d.\n", tokenCount);
                fclose(fp);
                return 1;
            }

            // 2) parse
            char  *tok, *save = NULL;
            int    n = 0;
            char   name[MAX_VAR_NAME];
            double fval = 0.0;
            for (tok = strtok_r(value, ",", &save);
                 tok;
                 tok = strtok_r(NULL, ",", &save), ++n)
            {
                 trim(tok);
                 if (n == 0) {
                     strncpy(name, tok, MAX_VAR_NAME-1);
                     name[MAX_VAR_NAME-1] = '\0';
                 } else if (n == 1) {
                     fval = atof(tok);
                 }
            }

            // 3) lookup & assign
            VariableBoundary *vb = findVariableBoundary(name, params);
            if (!vb) {
                fprintf(stderr, "Error: Fill_Constant for undeclared '%s'.\n", name);
                fclose(fp);
                return 1;
            }
            vb->fillType  = FILL_CONSTANT;
            vb->fillValue = fval;
            found_fill_constant = 1;

         
         } else if (strcasecmp(key, "RESPAWN") == 0) {
                params->RESPAWN = atoi(value);
                found_RESPAWN = 1;
            } else if (strcmp(key, "restart_time") == 0) {
                params->restart_time = atoi(value);
                found_restarttime = 1;
            } else if (strcasecmp(key, "WRITE_TO_CSV") == 0) {
                params->WRITE_TO_CSV = atoi(value);
                found_write_to_csv = 1;
            } else if (strcasecmp(key, "WRITE_TO_VTK") == 0) {
                params->WRITE_TO_VTK = atoi(value);
                found_write_to_vtk = 1;
            } else {
               // Unrecognized key: warn and skip.
               fprintf(stderr, "Warning: Unrecognized parameter '%s' in input file; skipping.\n", key);
           }
    }
    
    fclose(fp);

    if(params->DIM == 2) {
        params->Num_Z = 1;
        params->dz    = 1.0;
    } else if (params->DIM == 3) {
        /* require both Num_Z & dz */
        if (params->Num_Z <= 1) {
            fprintf(stderr, "Error: 'Num_Z' must be >1 for 3D.\n"); return 1;
        }
        if (params->dz <= 0.0) {
            fprintf(stderr, "Error: 'dz' must be >0 for 3D.\n"); return 1;
        }
    }
    
    /* Check for missing parameters and print appropriate messages */
    int error = 0;
    if (!found_DIM)            { fprintf(stderr, "Error: Parameter 'DIM' is missing.\n"); error = 1; }
    if (!found_Num_X)          { fprintf(stderr, "Error: Parameter 'Num_X' is missing.\n"); error = 1; }
    if (!found_Num_Y)          { fprintf(stderr, "Error: Parameter 'Num_Y' is missing.\n"); error = 1; }
    if (!found_dx)             { fprintf(stderr, "Error: Parameter 'dx' is missing.\n"); error = 1; }
    if (!found_dy)             { fprintf(stderr, "Error: Parameter 'dy' is missing.\n"); error = 1; }
    if (params->DIM == 3){
        if(!found_Num_Z)       { fprintf(stderr, "Error: Parameter 'Num_Z' is missing.\n"); error = 1;}
        if(!found_dz)          { fprintf(stderr, "Error: Parameter 'dz' is missing.\n"); error = 1;}
    }
    if (!found_dt)             { fprintf(stderr, "Error: Parameter 'dt' is missing.\n"); error = 1; }
    if (!found_total_steps)    { fprintf(stderr, "Error: Parameter 'total_steps' is missing.\n"); error = 1; }
    if (!found_timebreak)      { fprintf(stderr, "Error: Parameter 'timebreak' is missing.\n"); error = 1; }
    if (!found_epsilon)        { fprintf(stderr, "Error: Parameter 'epsilon' is missing.\n"); error = 1; }
    if (!found_tau)            { fprintf(stderr, "Error: Parameter 'tau' is missing.\n"); error = 1; }
    if (!found_delta)          { fprintf(stderr, "Error: Parameter 'delta' is missing.\n"); error = 1; }
    if (!found_j)              { fprintf(stderr, "Error: Parameter 'j' is missing.\n"); error = 1; }
    if (!found_theta_0)        { fprintf(stderr, "Error: Parameter 'theta_0' is missing.\n"); error = 1; }
    if (!found_alpha)          { fprintf(stderr, "Error: Parameter 'alpha' is missing.\n"); error = 1; }
    if (!found_gamma)          { fprintf(stderr, "Error: Parameter 'gamma' is missing.\n"); error = 1; }
    if (!found_a)              { fprintf(stderr, "Error: Parameter 'a' is missing.\n"); error = 1; }
    if (!found_K)              { fprintf(stderr, "Error: Parameter 'K' is missing.\n"); error = 1; }
    if (!found_T_e)            { fprintf(stderr, "Error: Parameter 'T_e' is missing.\n"); error = 1; }

    if (!found_boundary)       { fprintf(stderr, "Error: Parameter 'boundary' is missing.\n"); error = 1; }
    if ((!found_fill_cube) && (!found_fill_sphere) && (!found_fill_constant)) { fprintf(stderr, "Error: Filling info not provided.\n"); error = 1; }
    if ((found_RESPAWN) && (!found_restarttime)) { fprintf(stderr, "Error: Parameter 'restart_time' is missing.\n"); error = 1;}
    if ((!found_write_to_csv) && (!found_write_to_vtk)) { fprintf(stderr, "Error: File saving option not selected.\n"); error = 1;}

    if (error)
        return 1;
    
    return 0;  // success if all parameters were found
}



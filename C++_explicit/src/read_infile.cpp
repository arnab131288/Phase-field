#include "header.hpp"
#include <cstdio>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <errno.h>

/**
 * @brief Remove leading and trailing whitespace from a string in-place.
 *
 * @param str The input string to trim.
 */
void trim(char *str) {
    char *start = str;
    char *end;

    // Advance start past leading whitespace
    while (*start && std::isspace(static_cast<unsigned char>(*start))) {
        start++;
    }

    // Shift trimmed text back to the front
    if (start != str) {
        std::memmove(str, start, std::strlen(start) + 1);
    }

    // Trim trailing whitespace by null-terminating
    end = str + std::strlen(str) - 1;
    while (end >= str && std::isspace(static_cast<unsigned char>(*end))) {
        *end = '\0';
        end--;
    }
}

/**
 * @brief Read simulation parameters from a configuration file.
 *
 * Parses key=value; pairs and populates SimParams. Supports:
 *  - Grid dimensions and spacing (DIM, Num_X/Y/Z, dx/dy/dz)
 *  - Time stepping parameters (dt, total_steps, timebreak)
 *  - Material constants (epsilon, tau, delta, j, theta_0, alpha, gamma, a, K, T_e)
 *  - Boundary and fill specifications for each variable
 *  - Respawn and output options
 *
 * @param filename Path to the input file.
 * @param params   Pointer to SimParams to populate.
 * @return 0 on success, non-zero on error.
 */
int readParameters(const char *filename,  SimParams *params) {
    FILE *fp = std::fopen(filename, "r");
    if (!fp) {
        std::fprintf(stderr, "Error: Could not open file '%s': %s\n", filename, std::strerror(errno));
        return 1;
    }

    int found_DIM=0, found_Num_X=0, found_Num_Y=0, found_Num_Z=0;
    int found_dx=0, found_dy=0, found_dz=0, found_dt=0;
    int found_total_steps=0, found_timebreak=0;
    int found_epsilon=0, found_tau=0, found_delta=0, found_j=0;
    int found_theta_0=0, found_alpha=0, found_gamma=0;
    int found_a=0, found_K=0, found_T_e=0;
    int found_boundary=0, found_fill_cube=0, found_fill_sphere=0, found_fill_constant=0;
    int found_RESPAWN=0, found_restarttime=0;
    int found_write_to_csv=0, found_write_to_vtk=0;

    char line[256];
    while (std::fgets(line, sizeof(line), fp)) {
        if (line[0]=='#' || line[0]=='\n') continue;
        char key[64], value[64];
        if (!std::strchr(line,';') || std::sscanf(line, " %63[^=]=%63[^;];", key, value)!=2) {
            std::fprintf(stderr, "Warning: Could not parse line: %s", line);
            std::exit(EXIT_FAILURE);
        }
        trim(key); trim(value);

        if (strcasecmp(key,"DIM")==0)       { params->DIM=std::atoi(value); found_DIM=1; }
        else if (strcasecmp(key,"Num_X")==0) { params->Num_X=std::atoi(value); found_Num_X=1; }
        else if (strcasecmp(key,"Num_Y")==0) { params->Num_Y=std::atoi(value); found_Num_Y=1; }
        else if (strcasecmp(key,"Num_Z")==0) { params->Num_Z=std::atoi(value); found_Num_Z=1; }
        else if (strcasecmp(key,"dx")==0)    { params->dx=std::atof(value); found_dx=1; }
        else if (strcasecmp(key,"dy")==0)    { params->dy=std::atof(value); found_dy=1; }
        else if (strcasecmp(key,"dz")==0)    { params->dz=std::atof(value); found_dz=1; }
        else if (strcasecmp(key,"dt")==0)    { params->dt=std::atof(value); found_dt=1; }
        else if (strcasecmp(key,"total_steps")==0) { params->total_timesteps=std::atoi(value); found_total_steps=1; }
        else if (strcasecmp(key,"timebreak")==0)   { params->timebreak=std::atoi(value); found_timebreak=1; }
        else if (strcasecmp(key,"epsilon")==0) { params->epsilon=std::atof(value); found_epsilon=1; }
        else if (strcasecmp(key,"tau")==0)     { params->tau=std::atof(value); found_tau=1; }
        else if (strcasecmp(key,"delta")==0)   { params->delta=std::atof(value); found_delta=1; }
        else if (strcasecmp(key,"j")==0)       { params->j=std::atof(value); found_j=1; }
        else if (strcasecmp(key,"theta_0")==0){ params->theta_0=std::atof(value); found_theta_0=1; }
        else if (strcasecmp(key,"alpha")==0)  { params->alpha=std::atof(value); found_alpha=1; }
        else if (strcasecmp(key,"gamma")==0)  { params->gamma=std::atof(value); found_gamma=1; }
        else if (strcasecmp(key,"a")==0)      { params->a=std::atof(value); found_a=1; }
        else if (strcasecmp(key,"K")==0)      { params->K=std::atof(value); found_K=1; }
        else if (strcasecmp(key,"T_e")==0)    { params->T_e=std::atof(value); found_T_e=1; }
        else if (strcasecmp(key, "boundary") == 0) {
            // Expected format: variable,top,bottom,left,right[,front,back]
            char *token;
            int tokenCount = 0;
            char varName[MAX_VAR_NAME] = "";

            // Temporary FaceBoundary to hold the parsed boundary conditions.
            FaceBoundary tempBC;
            // Initialize to undefined
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
                    BoundaryType bcType;
                    if (strcasecmp(token, "NOFLUX") == 0) {
                        bcType = BOUNDARY_NOFLUX;
                    } else if (strcasecmp(token, "PERIODIC") == 0) {
                        bcType = BOUNDARY_PERIODIC;
                    } else if (strcasecmp(token, "UNDEFINED") == 0) {
                        bcType = BOUNDARY_UNDEFINED;
                    } else {
                        fprintf(stderr, "Error: unknown boundary condition '%s' for token %d.", token, tokenCount);
                        fclose(fp);
                        return 1;
                    }

                    switch (tokenCount) {
                        case 2: tempBC.top = bcType; break;
                        case 3: tempBC.bottom = bcType; break;
                        case 4: tempBC.left = bcType; break;
                        case 5: tempBC.right = bcType; break;
                        case 6: if (params->DIM == 3) tempBC.front = bcType; break;
                        case 7: if (params->DIM == 3) tempBC.back = bcType; break;
                        default:
                            fprintf(stderr, "Warning: Extra token '%s' in boundary specification for variable '%s'; ignoring.", token, varName);
                    }
                }
                token = strtok(NULL, ",");
            }

            if ((params->DIM == 2 && tokenCount != 5) || (params->DIM == 3 && tokenCount != 7)) {
                fprintf(stderr, "Error: Expected %d tokens for 'boundary' but got %d.",
                        (params->DIM == 2) ? 5 : 7, tokenCount);
                fclose(fp);
                return 1;
            }

            // Now update the SimParams structure.
            int found = 0;
            for (int i = 0; i < params->numVariables; i++) {
                if (strcasecmp(params->variables[i].varName, varName) == 0) {
                    params->variables[i].bc = tempBC;
                    found = 1;
                    break;
                }
            }
            if (!found) {
                if (params->numVariables >= MAX_VARIABLES) {
                    fprintf(stderr, "Error: Exceeded maximum number of variables (%d).", MAX_VARIABLES);
                    fclose(fp);
                    return 1;
                }
                strncpy(params->variables[params->numVariables].varName, varName, MAX_VAR_NAME - 1);
                params->variables[params->numVariables].varName[MAX_VAR_NAME - 1] = '\0';
                params->variables[params->numVariables].bc = tempBC;
                params->numVariables++;
            }
            found_boundary = 1;
        }
        else if (strcasecmp(key,"Fill_Cube")==0) {
            int n=0; char *save=NULL; char name[MAX_VAR_NAME]; double v=0; int x0=0,y0=0,z0=0,x1=0,y1=0,z1=0;
            for (char *t=strtok_r(value,",",&save); t; t=strtok_r(NULL,",",&save),n++) {
                trim(t);
                switch(n) {
                    case 0: strncpy(name,t,MAX_VAR_NAME-1); break;
                    case 1: v   = atof(t); break;
                    case 2: x0  = atoi(t); break;
                    case 3: x1  = atoi(t); break;
                    case 4: y0  = atoi(t); break;
                    case 5: y1  = atoi(t); break;
                    case 6: if (params->DIM==3) z0=atoi(t); break;
                    case 7: if (params->DIM==3) z1=atoi(t); break;
                }
            }
            VariableBoundary *vb = findVariableBoundary(name,params);
            vb->fillType = FILL_CUBE; vb->fillValue = v;
            vb->fill.c = {x0,y0,z0,x1,y1,z1}; found_fill_cube=1;
        }
        else if (strcasecmp(key,"Fill_Sphere")==0) {
            int n=0; char *save=NULL; char name[MAX_VAR_NAME]; double v=0; int r=0,cx=0,cy=0,cz=0;
            for (char *t=strtok_r(value,",",&save); t; t=strtok_r(NULL,",",&save),n++) {
                trim(t);
                switch(n) {
                    case 0: strncpy(name,t,MAX_VAR_NAME-1); break;
                    case 1: v      = atof(t); break;
                    case 2: r      = atoi(t); break;
                    case 3: cx     = atoi(t); break;
                    case 4: cy     = atoi(t); break;
                    case 5: if(params->DIM==3) cz=atoi(t); break;
                }
            }
            VariableBoundary *vb = findVariableBoundary(name,params);
            vb->fillType = FILL_SPHERE; vb->fillValue = v;
            vb->fill.s = {cx,cy,cz,r}; found_fill_sphere=1;
        }
        else if (strcasecmp(key,"Fill_Constant")==0) {
            int n=0; char *save=NULL; char name[MAX_VAR_NAME]; double v=0;
            for (char *t=strtok_r(value,",",&save); t; t=strtok_r(NULL,",",&save),n++) {
                trim(t);
                if(n==0) strncpy(name,t,MAX_VAR_NAME-1);
                else if(n==1) v=atof(t);
            }
            VariableBoundary *vb = findVariableBoundary(name,params);
            vb->fillType = FILL_CONSTANT; vb->fillValue = v; found_fill_constant=1;
        }
        else if (strcasecmp(key,"RESPAWN")==0)      { params->RESPAWN=atoi(value); found_RESPAWN=1; }
        else if (strcasecmp(key,"restart_time")==0){ params->restart_time=atoi(value); found_restarttime=1; }
        else if (strcasecmp(key,"WRITE_TO_CSV")==0){ params->WRITE_TO_CSV=atoi(value); found_write_to_csv=1; }
        else if (strcasecmp(key,"WRITE_TO_VTK")==0){ params->WRITE_TO_VTK=atoi(value); found_write_to_vtk=1; }
        else { std::fprintf(stderr,"Warning: Unrecognized key '%s'\n",key); }
    }
    std::fclose(fp);

    if (params->DIM==2) { params->Num_Z=1; params->dz=1.0; }
    else if (params->DIM==3 && (params->Num_Z<=1||params->dz<=0)) {
        std::fprintf(stderr,"Error: Invalid 3D parameters\n"); return 1;
    }

    int error=0;
    if(!found_DIM){fprintf(stderr,"Error: DIM missing.\n"); error=1;}    
    if(!found_Num_X){fprintf(stderr,"Error: Num_X missing.\n"); error=1;}    
    if(!found_Num_Y){fprintf(stderr,"Error: Num_Y missing.\n"); error=1;}    
    if(params->DIM==3 && !found_Num_Z){fprintf(stderr,"Error: Num_Z missing.\n"); error=1;}    
    if(!found_dx){fprintf(stderr,"Error: dx missing.\n"); error=1;}    
    if(!found_dy){fprintf(stderr,"Error: dy missing.\n"); error=1;}    
    if(params->DIM==3 && !found_dz){fprintf(stderr,"Error: dz missing.\n"); error=1;}    
    if(!found_dt){fprintf(stderr,"Error: dt missing.\n"); error=1;}    
    if(!found_total_steps){fprintf(stderr,"Error: total_steps missing.\n"); error=1;}    
    if(!found_timebreak){fprintf(stderr,"Error: timebreak missing.\n"); error=1;}    
    if(!found_epsilon){fprintf(stderr,"Error: epsilon missing.\n"); error=1;}    
    if(!found_tau){fprintf(stderr,"Error: tau missing.\n"); error=1;}    
    if(!found_delta){fprintf(stderr,"Error: delta missing.\n"); error=1;}    
    if(!found_j){fprintf(stderr,"Error: j missing.\n"); error=1;}    
    if(!found_theta_0){fprintf(stderr,"Error: theta_0 missing.\n"); error=1;}    
    if(!found_alpha){fprintf(stderr,"Error: alpha missing.\n"); error=1;}    
    if(!found_gamma){fprintf(stderr,"Error: gamma missing.\n"); error=1;}    
    if(!found_a){fprintf(stderr,"Error: a missing.\n"); error=1;}    
    if(!found_K){fprintf(stderr,"Error: K missing.\n"); error=1;}    
    if(!found_T_e){fprintf(stderr,"Error: T_e missing.\n"); error=1;}    
    if(!found_boundary){fprintf(stderr,"Error: boundary missing.\n"); error=1;}    
    if(!found_fill_cube&&!found_fill_sphere&&!found_fill_constant){fprintf(stderr,"Error: fill missing.\n"); error=1;}    
    if(found_RESPAWN&&!found_restarttime){fprintf(stderr,"Error: restart_time missing.\n"); error=1;}    
    if(!found_write_to_csv&&!found_write_to_vtk){fprintf(stderr,"Error: output option missing.\n"); error=1;}    

    return error?1:0;
}

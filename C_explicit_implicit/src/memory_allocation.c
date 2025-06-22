#include "header.h"

double *allocate_vector(int Num_X, int Num_Y, int Num_Z) {
    double *arr = (double *)malloc((Num_X * Num_Y * Num_Z) * sizeof(double));
    if(!arr) {
       fprintf(stderr, "Error: Could not allocate memory for 1D array.\n");
       exit(EXIT_FAILURE); 
    }
    return arr;
}


void free_vector(double *arr) {
    free(arr);
}

// Helper to allocate a contiguous 3D array of size NX*NY*NZ
static double *alloc3(int NX, int NY, int NZ) {
    double *ptr = allocate_vector(NX, NY, NZ);
    if (!ptr) {
        fprintf(stderr, "Error: allocate_vector failed for size %d x %d x %d\n", NX, NY, NZ);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void allocateFieldBuffers(const SimParams *params, FieldBuffers *fb) {
  int NX=params->Num_X, NY=params->Num_Y, NZ=params->Num_Z;
  fb->phi_new      = alloc3(NX,NY,NZ);
  fb->temp_new     = alloc3(NX,NY,NZ);
  fb->dphi_dt      = alloc3(NX,NY,NZ);
  fb->dfdphi       = alloc3(NX,NY,NZ);
  fb->ac           = alloc3(NX,NY,NZ);
  fb->ac_right     = alloc3(NX,NY,NZ);
  fb->ac_left      = alloc3(NX,NY,NZ);
  fb->ac_top       = alloc3(NX,NY,NZ);
  fb->ac_bottom    = alloc3(NX,NY,NZ);
  fb->ac_p         = alloc3(NX,NY,NZ);
  fb->ac_p_right   = alloc3(NX,NY,NZ);
  fb->ac_p_left    = alloc3(NX,NY,NZ);
  fb->ac_p_top     = alloc3(NX,NY,NZ);
  fb->ac_p_bottom  = alloc3(NX,NY,NZ); 
  fb->DERX_c       = alloc3(NX,NY,NZ);
  fb->DERY_c       = alloc3(NX,NY,NZ);
  fb->DERX_right   = alloc3(NX,NY,NZ);
  fb->DERX_left    = alloc3(NX,NY,NZ);
  fb->DERY_top     = alloc3(NX,NY,NZ);
  fb->DERY_bottom  = alloc3(NX,NY,NZ);
  fb->DERY_right   = alloc3(NX,NY,NZ);
  fb->DERY_left    = alloc3(NX,NY,NZ);
  fb->DERX_top     = alloc3(NX,NY,NZ);
  fb->DERX_bottom  = alloc3(NX,NY,NZ);
  fb->temp_buf1    = alloc3(NX,NY,NZ);
  fb->temp_buf2    = alloc3(NX,NY,NZ);
}

void freeFieldBuffers(FieldBuffers *fb) {
  free_vector(fb->phi_new);
  free_vector(fb->temp_new);
  free_vector(fb->dphi_dt);
  free_vector(fb->dfdphi);
  free_vector(fb->ac);
  free_vector(fb->ac_right);
  free_vector(fb->ac_left);
  free_vector(fb->ac_top);
  free_vector(fb->ac_bottom);
  free_vector(fb->ac_p);
  free_vector(fb->ac_p_right);
  free_vector(fb->ac_p_left);
  free_vector(fb->ac_p_top);
  free_vector(fb->ac_p_bottom);
  free_vector(fb->DERX_c);
  free_vector(fb->DERY_c);
  free_vector(fb->DERX_right);
  free_vector(fb->DERX_left);
  free_vector(fb->DERY_top);
  free_vector(fb->DERY_bottom);
  free_vector(fb->DERY_right);
  free_vector(fb->DERY_left);
  free_vector(fb->DERX_top);
  free_vector(fb->DERX_bottom);
  free_vector(fb->temp_buf1);
  free_vector(fb->temp_buf2);
}




// Define and initialize global storage for variable data.
VariableData globalVars[MAX_VARIABLES];
int numGlobalVars = 0;


void setupVariables(SimParams *params) {
    int totalElements;
    double *data;
    for (int i = 0; i < params->numVariables; i++) {
        char *varName = params->variables[i].varName;
        totalElements = params->Num_X * params->Num_Y * params->Num_Z;
        data = allocate_vector(params->Num_X, params->Num_Y, params->Num_Z);
        addVariableData(varName, data, totalElements);
    }
}

/*
 * Registers a variable and its data array into the global mapping.
 * varName: name of the variable
 * array: pointer to the data array
 * dataSize: the number of elements in the array (optional, can be used later)
 */
void addVariableData(const char* varName, double *array, size_t dataSize) {
    if (numGlobalVars >= MAX_VARIABLES) {
        fprintf(stderr, "Error: Exceeded maximum number of variables (%d).\n", MAX_VARIABLES);
        return;
    }
    // Copy the variable name and store the pointer.
    strncpy(globalVars[numGlobalVars].varName, varName, MAX_VAR_NAME - 1);
    globalVars[numGlobalVars].varName[MAX_VAR_NAME - 1] = '\0';
    globalVars[numGlobalVars].dataArray = array;
    globalVars[numGlobalVars].dataSize = dataSize;
    numGlobalVars++;
}

/*
 * Returns the data array pointer associated with the given variable name.
 * If the variable is not found, returns NULL.
 */
double* getDataArray(const char* varName) {
    for (int i = 0; i < numGlobalVars; i++) {
        if (strcmp(globalVars[i].varName, varName) == 0) {
            return globalVars[i].dataArray;
        }
    }
    fprintf(stderr, "Warning: Variable '%s' not found.\n", varName);
    return NULL;
}

VariableBoundary* findVariableBoundary(const char *name, SimParams *params) {
    for (int i = 0; i < params->numVariables; i++) {
        if (strcmp(params->variables[i].varName, name) == 0) {
            return &(params->variables[i]);
        }
    }
    return NULL;
}


void freeGlobalVariableArrays(void) {
    for (int i = 0; i < numGlobalVars; i++) {
        free(globalVars[i].dataArray);
        globalVars[i].dataArray = NULL;  // Optional: clear pointer after free.
    }
    numGlobalVars = 0;  // Optionally reset the count.
}

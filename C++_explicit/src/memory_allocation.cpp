#include "header.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>

/*
 * memory_allocation.cpp
 *
 * Provides routines for dynamic allocation, deallocation, and global management
 * of simulation data arrays:
 *  - allocate_vector: allocate a flat 1D array for NX*NY*NZ elements via malloc
 *  - free_vector: free memory allocated by allocate_vector
 *  - alloc3: helper to allocate contiguous 3D data via allocate_vector
 *  - allocateFieldBuffers: allocate all intermediate FieldBuffers arrays
 *  - freeFieldBuffers: release all memory in FieldBuffers
 *  - setupVariables: allocate and register each simulation variable
 *  - addVariableData: register a variable and its data array globally
 *  - getDataArray: retrieve a data array by variable name
 *  - findVariableBoundary: find boundary settings by variable name
 *  - freeGlobalVariableArrays: free all registered arrays
 */

/**
 * @brief Allocate a flat vector of doubles for a 3D grid using malloc.
 * Exits on failure.
 */
double* allocate_vector(int Num_X, int Num_Y, int Num_Z) {
    size_t total = static_cast<size_t>(Num_X) * Num_Y * Num_Z;
    double *arr = static_cast<double*>(std::malloc(total * sizeof(double)));
    if (!arr) {
        std::fprintf(stderr, "Error: Could not allocate memory for %zu elements.\n", total);
        std::exit(EXIT_FAILURE);
    }
    return arr;
}

/**
 * @brief Free a vector allocated by allocate_vector.
 */
void free_vector(double *arr) {
    std::free(arr);
}

/**
 * @brief Helper to allocate a contiguous 3D block of data.
 */
double* alloc3(int NX, int NY, int NZ) {
    double* ptr = allocate_vector(NX, NY, NZ);
    return ptr;
}

/**
 * @brief Allocate all intermediate buffers in a FieldBuffers struct.
 */
void allocateFieldBuffers(const SimParams *params, FieldBuffers *fb) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    fb->phi_new      = alloc3(NX, NY, NZ);
    fb->temp_new     = alloc3(NX, NY, NZ);
    fb->dphi_dt      = alloc3(NX, NY, NZ);
    fb->dfdphi       = alloc3(NX, NY, NZ);
    fb->ac           = alloc3(NX, NY, NZ);
    fb->ac_right     = alloc3(NX, NY, NZ);
    fb->ac_left      = alloc3(NX, NY, NZ);
    fb->ac_top       = alloc3(NX, NY, NZ);
    fb->ac_bottom    = alloc3(NX, NY, NZ);
    fb->ac_p         = alloc3(NX, NY, NZ);
    fb->ac_p_right   = alloc3(NX, NY, NZ);
    fb->ac_p_left    = alloc3(NX, NY, NZ);
    fb->ac_p_top     = alloc3(NX, NY, NZ);
    fb->ac_p_bottom  = alloc3(NX, NY, NZ);
    fb->DERX_c       = alloc3(NX, NY, NZ);
    fb->DERY_c       = alloc3(NX, NY, NZ);
    fb->DERX_right   = alloc3(NX, NY, NZ);
    fb->DERX_left    = alloc3(NX, NY, NZ);
    fb->DERY_top     = alloc3(NX, NY, NZ);
    fb->DERY_bottom  = alloc3(NX, NY, NZ);
    fb->DERY_right   = alloc3(NX, NY, NZ);
    fb->DERY_left    = alloc3(NX, NY, NZ);
    fb->DERX_top     = alloc3(NX, NY, NZ);
    fb->DERX_bottom  = alloc3(NX, NY, NZ);
}

/**
 * @brief Free all intermediate buffers in a FieldBuffers struct.
 */
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
}

// Define and initialize global storage for variable data.
VariableData globalVars[MAX_VARIABLES];
int numGlobalVars = 0;

/**
 * @brief Allocate and register each variable array based on SimParams.
 */
void setupVariables(const SimParams *params) {
    size_t totalElements = static_cast<size_t>(params->Num_X) * params->Num_Y * params->Num_Z;
    for (int i = 0; i < params->numVariables; ++i) {
        const char* name = params->variables[i].varName;
        double* data = allocate_vector(params->Num_X, params->Num_Y, params->Num_Z);
        addVariableData(name, data, totalElements);
    }
}

/**
 * @brief Register a variable and its data array.
 */
void addVariableData(const char* varName, double *array, size_t dataSize) {
    if (numGlobalVars >= MAX_VARIABLES) {
        std::fprintf(stderr, "Error: Exceeded maximum number of variables (%d).\n", MAX_VARIABLES);
        return;
    }
    std::strncpy(globalVars[numGlobalVars].varName, varName, MAX_VAR_NAME - 1);
    globalVars[numGlobalVars].varName[MAX_VAR_NAME - 1] = '\0';
    globalVars[numGlobalVars].dataArray = array;
    globalVars[numGlobalVars].dataSize = dataSize;
    ++numGlobalVars;
}

/**
 * @brief Retrieve the data array pointer for a given variable name.
 */
double* getDataArray(const char* varName) {
    for (int i = 0; i < numGlobalVars; ++i) {
        if (std::strcmp(globalVars[i].varName, varName) == 0) {
            return globalVars[i].dataArray;
        }
    }
    std::fprintf(stderr, "Warning: Variable '%s' not found.\n", varName);
    return nullptr;
}

/**
 * @brief Find the VariableBoundary for a given variable name in SimParams.
 */
VariableBoundary* findVariableBoundary(const char *name, SimParams *params) {
    for (int i = 0; i < params->numVariables; ++i) {
        if (std::strcmp(params->variables[i].varName, name) == 0) {
            return &params->variables[i];
        }
    }
    return nullptr;
}

/**
 * @brief Free all globally registered variable arrays.
 */
void freeGlobalVariableArrays(void) {
    for (int i = 0; i < numGlobalVars; ++i) {
        std::free(globalVars[i].dataArray);
        globalVars[i].dataArray = nullptr;
    }
    numGlobalVars = 0;
}

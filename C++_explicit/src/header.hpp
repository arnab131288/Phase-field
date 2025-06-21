#ifndef HEADER_HPP
#define HEADER_HPP

/*
 * header.hpp
 *
 * Declares the core data structures, constants, and function prototypes
 * for a multidimensional simulation framework. Defines:
 *  - Simulation parameters (SimParams)
 *  - Variable boundary and filling types (BoundaryType, FillType, VariableBoundary)
 *  - Field buffer containers for intermediate computations (FieldBuffers)
 *  - Variable data mapping (VariableData, globalVars)
 *  - Allocation/deallocation routines for buffers
 *  - I/O functions for parameters, VTK/CSV input and output
 *  - Simulation routines: filling shapes, updating fields, computing derivatives
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <cerrno>
#include <unistd.h>

//-----------------------------------------------------------------------------
// Constants and limits
//----------------------------------------------------------------------------- 
static constexpr int MAX_VAR_NAME    = 32;   // Maximum length for variable names.
static constexpr int MAX_VARIABLES   = 10;   // Maximum number of variables supported.
static constexpr int MAX_DIM         = 3;    // Maximum spatial dimensions     

//-----------------------------------------------------------------------------
// Enum to represent boundary & filling type.
//----------------------------------------------------------------------------- 
enum BoundaryType {
    BOUNDARY_UNDEFINED,
    BOUNDARY_NOFLUX,
    BOUNDARY_PERIODIC
    // Additional boundary conditions can be added here.
};

enum FillType { 
    FILL_NONE, 
    FILL_CUBE, 
    FILL_SPHERE,
    FILL_CONSTANT 
    // Additional filling types can be added here.
};

//-----------------------------------------------------------------------------
// Boundary definitions per face
//----------------------------------------------------------------------------- 
struct FaceBoundary {
    BoundaryType top;
    BoundaryType bottom;
    BoundaryType left;
    BoundaryType right;
    BoundaryType front;
    BoundaryType back;
};

//-----------------------------------------------------------------------------
// Geometric shapes for filling domains
//----------------------------------------------------------------------------- 
struct cube {
    int x_start, y_start, z_start;
    int x_end,   y_end,   z_end;
};

struct sphere {
    int x_center, y_center, z_center;
    int radius;
};

//-----------------------------------------------------------------------------
// Variable boundary and data structures
//----------------------------------------------------------------------------- 
struct VariableBoundary {
    char varName[MAX_VAR_NAME];
    FaceBoundary bc;
    FillType    fillType;
    union {
        cube   c;
        sphere s;
    } fill;
    double fillValue;
};

struct VariableData {
    char   varName[MAX_VAR_NAME];
    double *dataArray;  // Pointer to the variable's data
    size_t  dataSize;   // Number of elements in the array
};

//-----------------------------------------------------------------------------
// Struct to hold simulation parameters.
//----------------------------------------------------------------------------- 
struct SimParams {
    int    DIM;
    int    Num_X;
    int    Num_Y;
    int    Num_Z;
    double dx, dy, dz;
    double dt;
    int    total_timesteps;
    int    timebreak;
    double epsilon, tau, delta;
    int    j;
    double theta_0, alpha, gamma, a, K, T_e;
    int    numVariables;
    VariableBoundary variables[MAX_VARIABLES];

    // Respawn parameters
    int RESPAWN;
    int restart_time;

    // File writing options
    int WRITE_TO_CSV;
    int WRITE_TO_VTK;
};

//-----------------------------------------------------------------------------
// Field buffers for intermediate computations
//----------------------------------------------------------------------------- 
struct FieldBuffers {
    double *phi_new, *temp_new;
    double *dphi_dt, *dfdphi;
    double *ac, *ac_right, *ac_left, *ac_top, *ac_bottom;
    double *ac_p, *ac_p_right, *ac_p_left, *ac_p_top, *ac_p_bottom;
    double *DERX_c, *DERY_c;
    double *DERX_right, *DERX_left, *DERY_top, *DERY_bottom;
    double *DERY_right, *DERY_left, *DERX_top, *DERX_bottom;
};

//-----------------------------------------------------------------------------
// Globals for external variable data mapping
//----------------------------------------------------------------------------- 
extern VariableData globalVars[MAX_VARIABLES];
extern int          numGlobalVars;

//-----------------------------------------------------------------------------
// Allocation / deallocation routines
//----------------------------------------------------------------------------- 
void allocateFieldBuffers(const SimParams *params, FieldBuffers *fb);
void freeFieldBuffers(FieldBuffers *fb);

//-----------------------------------------------------------------------------
// Function prototypes for file I/O.
//----------------------------------------------------------------------------- 
int    readParameters(const char *filename,  SimParams *params);
void   writeParameters(const char *outfile, const SimParams *params);
void   read_input_vtk(const char *filename, double *arr, const SimParams *params, int strides[]);
void   read_input_csv(const char *filename, double *arr, const SimParams *params, int strides[]);
void   write_output_vtk(const char *filename, double *arr, const SimParams *params, int strides[]);
void   write_output_csv(const char *filename, double *arr, const SimParams *params, int strides[]);
void   trim(char *str);

//-----------------------------------------------------------------------------
// Variable management routines
//----------------------------------------------------------------------------- 
void   setupVariables(const SimParams *params);
double *alloc3(int NX, int NY, int NZ);
void   freeGlobalVariableArrays(void);
void   addVariableData(const char* varName, double *array, size_t dataSize);
double* getDataArray(const char* varName);    
VariableBoundary* findVariableBoundary(const char *name, SimParams *params);

//-----------------------------------------------------------------------------
// Function prototype for applying boundary conditions.
//----------------------------------------------------------------------------- 
void applyBoundaryConditions(double *arr, const SimParams *params, int strides[], FaceBoundary bc);

//-----------------------------------------------------------------------------
// Function prototypes for simulation routines.
//----------------------------------------------------------------------------- 
void   FillCube(double *arr, const VariableBoundary *vb, const SimParams *params, int strides[]);
void   FillSphere(double *arr, const VariableBoundary *vb, const SimParams *params, int strides[]);
void   FillConstant(double *arr, const VariableBoundary *vb, const SimParams *params, int strides[]);
void   updateTemp(double *temp, FieldBuffers *fb, const SimParams *params, int strides[], double r2[]);
double computeLaplacian(double *arr, int index, int strides[], double r2[], int dim);
void   computedfdphi(double *phi, double *dfdphi, double *temp, const SimParams *params, int strides[]);
void   updatePhi(double *phi, FieldBuffers *fb, const SimParams *params, double r[], int strides[]);
void   computeGradientPhi(double *phi, FieldBuffers *fb, const SimParams *params, double r[], int strides[]); 
void   computeAnisotropy(FieldBuffers *fb, const SimParams *params, int strides[]);     
void   copyInterior(double *dst, double *src, const SimParams *params, int strides[]);

#endif // HEADER_HPP

#ifndef HEADER_H
#define HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>

//-----------------------------------------------------------------------------
// Constants and limits
//-----------------------------------------------------------------------------
#define MAX_VAR_NAME 32      // Maximum length for variable names.
#define MAX_VARIABLES 10     // Maximum number of variables supported (can be increased dynamically as needed).
#define MAX_DIM  3           // Maximum spatial dimensions     

/*--------------------------------------------------------------
 * Enum to represent boundary & filling type.
 *--------------------------------------------------------------*/
typedef enum {
    BOUNDARY_UNDEFINED,
    BOUNDARY_NOFLUX,
    BOUNDARY_PERIODIC
    // Additional boundary conditions can be added here.
} BoundaryType;

typedef enum { 
    FILL_NONE, 
    FILL_CUBE, 
    FILL_SPHERE,
    FILL_CONSTANT 
    //Additional filling type can be added here.
} FillType;

//-----------------------------------------------------------------------------
// Boundary definitions per face
//-----------------------------------------------------------------------------
typedef struct {
    BoundaryType top;
    BoundaryType bottom;
    BoundaryType left;
    BoundaryType right;
    BoundaryType front;
    BoundaryType back;
} FaceBoundary;

//-----------------------------------------------------------------------------
// Geometric shapes for filling domains
//-----------------------------------------------------------------------------
typedef struct {
  int x_start, y_start, z_start;
  int x_end,   y_end,   z_end;
} cube;

typedef struct {
  int x_center, y_center, z_center;
  int radius;
} sphere;

//-----------------------------------------------------------------------------
// Variable boundary and data structures
//-----------------------------------------------------------------------------


typedef struct {
    char varName[MAX_VAR_NAME];
    FaceBoundary bc;
    FillType    fillType;
    union {
      cube   c;
      sphere s;
  } fill;
  double fillValue;
} VariableBoundary;

typedef struct {
    char varName[MAX_VAR_NAME];
    double *dataArray;  // Pointer to the variable's data
    size_t dataSize;    // Optionally, store the size (number of elements) in the array
} VariableData;

/*--------------------------------------------------------------
 * Struct to hold simulation parameters.
 *--------------------------------------------------------------*/
typedef struct SimParams{
    int DIM;
    int Num_X;
    int Num_Y;
    int Num_Z;
    double dx;
    double dy;
    double dz;
    double dt;
    int total_timesteps;
    int timebreak;
    double epsilon;
    double tau;
    double delta;
    int j;
    double theta_0;
    double alpha;
    double gamma;
    double a;
    double K;
    double T_e;
    int numVariables;                                
    VariableBoundary variables[MAX_VARIABLES];  

    /* Respawn parameters*/
    int RESPAWN;
    int restart_time;

    /*File writing options*/
    int WRITE_TO_CSV;
    int WRITE_TO_VTK;
} SimParams;

//-----------------------------------------------------------------------------
// Field buffers for intermediate computations
//-----------------------------------------------------------------------------

typedef struct FieldBuffers{
  double *phi_new;
  double *temp_new;
  double *dphi_dt;
  double *dfdphi;
  double *ac;
  double *ac_right;
  double *ac_left;
  double *ac_top;
  double *ac_bottom;
  double *ac_p;
  double *ac_p_right;
  double *ac_p_left;
  double *ac_p_top;
  double *ac_p_bottom;
  double *DERX_c;
  double *DERY_c;
  double *DERX_right;
  double *DERX_left;
  double *DERY_top;
  double *DERY_bottom;
  double *DERY_right;
  double *DERY_left;
  double *DERX_top;
  double *DERX_bottom;
  double *temp_buf1;
  double *temp_buf2;
} FieldBuffers;



//-----------------------------------------------------------------------------
// Globals for external variable data mapping
//-----------------------------------------------------------------------------
extern VariableData globalVars[MAX_VARIABLES];
extern int numGlobalVars;


//-----------------------------------------------------------------------------
// Allocation / deallocation routines
//-----------------------------------------------------------------------------

void allocateFieldBuffers(const SimParams *params, FieldBuffers *fb);
void freeFieldBuffers(FieldBuffers *fb);



/*--------------------------------------------------------------
 * Function prototypes for file I/O.
 *--------------------------------------------------------------*/
int readParameters(const char *filename, SimParams *params);
void writeParameters(const char *outfile, const SimParams *params);
void read_input_vtk(const char *filename, double *arr, SimParams *params, int strides[]);
void read_input_csv(const char *filename, double *arr, SimParams *params, int strides[]);
void write_output_vtk(const char *filename, double *arr, SimParams *params, int strides[]);
void write_output_csv(const char *filename, double *arr, SimParams *params, int strides[]);
void trim(char *str);

/*--------------------------------------------------------------
 * Variable management routines
 *--------------------------------------------------------------*/
void setupVariables(SimParams *params);
static double *alloc3(int NX, int NY, int NZ);
void freeGlobalVariableArrays(void);
void addVariableData(const char* varName, double *array, size_t dataSize);
double* getDataArray(const char* varName);    
VariableBoundary* findVariableBoundary(const char *name, SimParams *params);

/*--------------------------------------------------------------
 * Function prototype for applying boundary conditions.
 *--------------------------------------------------------------*/
 void applyBoundaryConditions(double *arr, SimParams *params, int strides[], FaceBoundary bc);

/*--------------------------------------------------------------
 * Function prototypes for simulation routines.
 *--------------------------------------------------------------*/
void FillCube(double *arr, const VariableBoundary *vb, SimParams *params, int strides[]);
void FillSphere(double *arr, const VariableBoundary *vb, SimParams *params, int strides[]);
void FillConstant(double *arr, const VariableBoundary *vb, SimParams *params, int strides[]);
void updateTemp(double *temp, FieldBuffers *fb, SimParams *params, int strides[], double r2[]);
double computeLaplacian(double *arr, int index, int strides[], double r2[], int dim);
void computedfdphi(double *phi, double *dfdphi, double *temp, SimParams *params, int strides[]);
void updatePhi(double *phi, FieldBuffers *fb, SimParams *params, double r[], int strides[]);
void computeGradientPhi(double *phi, FieldBuffers *fb, SimParams *params, double r[], int strides[]); 
void computeAnisotropy(FieldBuffers *fb, SimParams *params, int strides[]);     
void copyInterior(double *dst, double *src, SimParams *params, int strides[]);

                 

#endif /* HEADER_H */


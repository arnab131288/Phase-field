#include "header.h"

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])

// discrete second derivative in X direction (no spacing factor)
static inline double lapX(double *field, int i, int j, int *strides) {
    return field[IDX(i+1,j,0)] - 2*field[IDX(i,j,0)] + field[IDX(i-1,j,0)];
}

// discrete second derivative in Y direction (no spacing factor)
static inline double lapY(double *field, int i, int j, int *strides) {
    return field[IDX(i,j+1,0)] - 2*field[IDX(i,j,0)] + field[IDX(i,j-1,0)];
}


// Solve tridiagonal system a*x_{k-1} + b*x_k + c*x_{k+1} = d_k of size M
static void thomas(int M, double *a, double *b, double *c, double *d, double *x) {
    double *c_prime = malloc(M * sizeof(double));
    double *d_prime = malloc(M * sizeof(double));
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for (int p = 1; p < M; p++) {
        double m = b[p] - a[p] * c_prime[p - 1];
        c_prime[p] = c[p] / m;
        d_prime[p] = (d[p] - a[p] * d_prime[p - 1]) / m;
    }
    x[M - 1] = d_prime[M - 1];
    for (int p = M - 2; p >= 0; p--) {
        x[p] = d_prime[p] - c_prime[p] * x[p + 1];
    }
    free(c_prime);
    free(d_prime);
}

// ADI step 1: implicit in X, explicit in Y using discrete Laplacians
static void adi_step1_X(int NX, int NY, double rx, double ry,
                        double *in, double *out, int strides[]) {
    int M = NX - 2;  // system size
    double *a = malloc(M * sizeof(double));
    double *b = malloc(M * sizeof(double));
    double *c = malloc(M * sizeof(double));
    double *d = malloc(M * sizeof(double));
    double *x = malloc(M * sizeof(double));

    for (int j = 1; j < NY - 1; j++) {
        // build tridiagonal for row j
        for (int i = 1; i < NX - 1; i++) {
            int k = i - 1;
            a[k] = -rx;
            b[k] = 1.0 + 2.0*rx;
            c[k] = -rx;
            d[k] = in[IDX(i,j,0)] + ry * lapY(in, i, j, strides);
        }
        // Neumann BC adjustments
        b[0]     -= a[0];
        b[M - 1] -= c[M - 1];
        thomas(M, a, b, c, d, x);
        // write back solution
        for (int i = 1; i < NX - 1; i++) {
            out[IDX(i,j,0)] = x[i - 1];
        }
    }
    free(a); free(b); free(c); free(d); free(x);
}

// ADI step 2: implicit in Y, explicit in X using discrete Laplacians
static void adi_step2_Y(int NX, int NY, double rx, double ry,
                        double *in, double *out, int strides[]) {
    int M = NY - 2;
    double *a = malloc(M * sizeof(double));
    double *b = malloc(M * sizeof(double));
    double *c = malloc(M * sizeof(double));
    double *d = malloc(M * sizeof(double));
    double *x = malloc(M * sizeof(double));

    for (int i = 1; i < NX - 1; i++) {
        // build tridiagonal for column i
        for (int j = 1; j < NY - 1; j++) {
            int p = j - 1;
            a[p] = -ry;
            b[p] = 1 + 2*ry;
            c[p] = -ry;
            d[p] = in[IDX(i,j,0)] + rx * lapX(in, i, j, strides);
        }
        // Neumann BC adjustments
        b[0]     -= a[0];
        b[M - 1] -= c[M - 1];
        thomas(M, a, b, c, d, x);
        // write back solution
        for (int j = 1; j < NY - 1; j++) {
            out[IDX(i,j,0)] = x[j - 1];
        }
    }
    free(a); free(b); free(c); free(d); free(x);
}

void updateTemp(double *T, FieldBuffers *fb, SimParams *params, int strides[], double r2[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    double dx = params->dx, dy = params->dy;
    double dt = params->dt;
    double K = params->K;

    double rx = 0.5 * K * dt * r2[0];
    double ry = 0.5 * K * dt * r2[1];

    double *T0 = T;
    double *T1 = fb->temp_buf1;
    double *T2 = fb->temp_buf2;
    double *Tn = fb->temp_new;
    
    // inherit phi BC on the dphi_dt buffer, so its ghost-cells are valid
    VariableBoundary *vb_phi = findVariableBoundary("phi", params);
    if (vb_phi) {
        applyBoundaryConditions(fb->dphi_dt, params, strides, vb_phi->bc);
    }

    
    for (int i = 1; i < NX - 1; i++) {
        for(int j = 1; j < NY - 1; j++) {
            int index = IDX(i, j, 0);
            T1[index] = T0[index] + 0.5 * dt * K * fb->dphi_dt[index];
        }
    }

    VariableBoundary *vb_temp = findVariableBoundary("temp", params);
    if (vb_temp) {
           applyBoundaryConditions(T1, params, strides, vb_temp->bc);
    }


    adi_step1_X(NX, NY, rx, ry, T1, T2, strides);

         
    if (vb_temp) {
        applyBoundaryConditions(T2, params, strides, vb_temp->bc);
    }


    for(int i = 1; i < NX - 1; i++) {
        for(int j = 1; j < NY - 1; j++) {
            int index = IDX(i,j,0);
            T2[index] += 0.5 * dt * K * fb->dphi_dt[index];
        }
    }
    //applyTemperatureBC(T2, params, strides);
    adi_step2_Y(NX, NY, rx, ry, T2, T1, strides);
    if (vb_temp) {
       applyBoundaryConditions(T1, params, strides, vb_temp->bc);
    }
    

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            int index = IDX(i, j, 0);
            Tn[index] = T1[index];
        }
    }
        
}

#undef IDX

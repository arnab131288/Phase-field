#include "header.h"
#include<math.h>

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])

void computedfdphi(double *phi, double *dfdphi, double *temp, SimParams *params, int strides[]) {
     int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
     int kstart = (dim == 3) ? 1 : 0;
     int kend   = (dim == 3) ? NZ - 1 : 1;

     for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            int index = IDX(i, j, 0);
            double m = ((params->alpha)/M_PI)*atan(params->gamma*(params->T_e - temp[index]));
            dfdphi[index] = phi[index] * (1.0 - phi[index]) * (phi[index] - 0.5 + m); 
        }
     }
}



#undef IDX

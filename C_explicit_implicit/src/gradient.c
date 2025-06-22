#include "header.h"

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])

void computeGradientPhi(double *phi, FieldBuffers *fb, SimParams *params, double r[], int strides[]) {
     int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
     double dx = params->dx, dy = params->dy, dz = params->dz;
     int kstart = (dim == 3) ? 1 : 0;
     int kend   = (dim == 3) ? NZ - 1 : 1;

     for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            int index = IDX(i, j, 0);

            /*Derivatives of phase-field via Forward difference*/
            fb->DERX_right[index]     = (phi[index + strides[0]] - phi[index]) * r[0];
			fb->DERX_left[index]      = (phi[index]              - phi[index - strides[0]]) * r[0];
			fb->DERY_top[index]       = (phi[index + strides[1]] - phi[index]) * r[1];
			fb->DERY_bottom[index]    = (phi[index]              - phi[index - strides[1]]) * r[1];

            /*Derivatives of phase-field via Central difference*/
			fb->DERX_c[index]         = 0.5 * (phi[index + strides[0]] - phi[index - strides[0]]) * r[0];
			fb->DERY_c[index]         = 0.5 * (phi[index + strides[1]] - phi[index - strides[1]]) * r[1];

            /*X derivative in the Y direction (top & bottom) and Y derivative in the X direction (right & left) using neighbours and next-neighbours*/
			fb->DERX_top[index]     = 0.25 * ((phi[index + strides[0]] + phi[index + strides[0] + strides[1]] + phi[index] + phi[index + strides[1]]) - (phi[index - strides[0]] + phi[index - strides[0] + strides[1]] + phi[index] + phi[index + strides[1]])) * r[0];
			fb->DERX_bottom[index]  = 0.25 * ((phi[index + strides[0]] + phi[index + strides[0] - strides[1]] + phi[index] + phi[index - strides[1]]) - (phi[index - strides[0]] + phi[index - strides[0] - strides[1]] + phi[index] + phi[index - strides[1]])) * r[0];
			fb->DERY_right[index]   = 0.25 * ((phi[index + strides[1]] + phi[index + strides[0] + strides[1]] + phi[index] + phi[index + strides[0]]) - (phi[index - strides[1]] + phi[index + strides[0] - strides[1]] + phi[index] + phi[index + strides[0]])) * r[1];
			fb->DERY_left[index]    = 0.25 * ((phi[index + strides[1]] + phi[index - strides[0] + strides[1]] + phi[index] + phi[index - strides[0]]) - (phi[index - strides[1]] + phi[index - strides[0] - strides[1]] + phi[index] + phi[index - strides[0]])) * r[1];
		}
	 }

}



#undef IDX

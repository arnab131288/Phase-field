#include "header.h"

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])


void updatePhi(double *phi, FieldBuffers *fb, SimParams *params, double r[], int strides[]) {
    
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;
    double dx = params->dx, dy = params->dy, dz = params->dz;
    double rj, lj, tj, bj, noise;
    int kstart = (dim == 3) ? 1 : 0;
    int kend   = (dim == 3) ? NZ - 1 : 1;

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            int index = IDX(i, j, 0);

            /*Phase-field fluxes on the right, left, top and bottom edges*/
            rj = fb->ac_right[index] * (fb->ac_right[index] * fb->DERX_right[index]  - fb->ac_p_right[index] * fb->DERY_right[index]);
	        lj = fb->ac_left[index]  * (fb->ac_left[index]  * fb->DERX_left[index]   - fb->ac_p_left[index]  * fb->DERY_left[index]);
	        tj = fb->ac_top[index]   * (fb->ac_top[index]   * fb->DERY_top[index]    + fb->ac_p_top[index]   * fb->DERX_top[index]);
	        bj = fb->ac_bottom[index]* (fb->ac_bottom[index]* fb->DERY_bottom[index] + fb->ac_p_bottom[index]* fb->DERX_bottom[index]);

            noise  = params->a*((double) rand()/RAND_MAX - 0.5); 
            noise *= phi[index] * (1.0 - phi[index]);

            fb->dphi_dt[index]  = (rj-lj) * r[0] + (tj-bj) * r[1];
            fb->dphi_dt[index]  += fb->dfdphi[index];
            fb->dphi_dt[index]  += noise;
            fb->dphi_dt[index]  /= params->tau; 

            fb->phi_new[index]  = phi[index] + params->dt*fb->dphi_dt[index];

        }
    }
     				
}

void copyInterior(double *dst, double *src, SimParams *params, int strides[]) {
    int NX = params->Num_X;
    int NY = params->Num_Y;
    int NZ = params->Num_Z;
    int dim = params->DIM;

        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                int index = IDX(i, j, 0);
                dst[index] = src[index];
            }
        }
    
}

#undef IDX

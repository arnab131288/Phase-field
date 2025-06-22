#include "header.h"
#include<math.h>

#define IDX(i,j,k) ((i)*strides[0] + (j)*strides[1] + (k)*strides[2])


 void computeAnisotropy(FieldBuffers *fb, SimParams *params, int strides[]) {
	 int NX = params->Num_X, NY = params->Num_Y, NZ = params->Num_Z, dim = params->DIM;
     double dx = params->dx, dy = params->dy, dz = params->dz, eps = params->epsilon;
	 double theta, theta_right, theta_left, theta_bottom, theta_top;
     int kstart = (dim == 3) ? 1 : 0;
     int kend   = (dim == 3) ? NZ - 1 : 1;

     for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            int index = IDX(i, j, 0);

			/*Anisotropy functions and derivative*/
			theta  		              = atan2(fb->DERY_c[index],fb->DERX_c[index]);
			fb->ac[index]             =  eps * (1.0 + params->delta * cos(params->j * (theta - params->theta_0)));
			fb->ac_p[index]           = -eps *(params->delta * params->j * sin(params->j * (theta - params->theta_0)));

			/*Normal to the interface*/
			theta_right           = atan2(fb->DERY_right[index],   fb->DERX_right[index]);
			theta_left            = atan2(fb->DERY_left[index],    fb->DERX_left[index]);
			theta_top             = atan2(fb->DERY_top[index],     fb->DERX_top[index]);
			theta_bottom          = atan2(fb->DERY_bottom[index],  fb->DERX_bottom[index]);

            /*Anisotropy function on neighbouring grid points*/
			fb->ac_right[index]		  =  eps * (1.0 + params->delta * cos(params->j * (theta_right-params->theta_0)));
			fb->ac_left[index]		  =  eps * (1.0 + params->delta * cos(params->j * (theta_left-params->theta_0)));
			fb->ac_top[index]    	  =  eps * (1.0 + params->delta * cos(params->j * (theta_top-params->theta_0)));
			fb->ac_bottom[index] 	  =  eps * (1.0 + params->delta * cos(params->j * (theta_bottom-params->theta_0)));

            /*Derivative of anisotropy function on neighbouring grid points*/
			fb->ac_p_right[index] 	  = -eps * (params->delta * params->j * sin(params->j * (theta_right-params->theta_0)));
			fb->ac_p_left[index]  	  = -eps * (params->delta * params->j * sin(params->j * (theta_left-params->theta_0)));
			fb->ac_p_top[index]   	  = -eps * (params->delta * params->j * sin(params->j * (theta_top-params->theta_0)));
			fb->ac_p_bottom[index] 	  = -eps * (params->delta * params->j * sin(params->j * (theta_bottom-params->theta_0)));

		}
	}
}



#undef IDX

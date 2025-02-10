#ifndef _COMPUTEMODELS_H
#define _COMPUTEMODELS_H

#include "mkl_types.h"

typedef struct
{
	// HOST_WORLD, GLOBAL DATA
	double cell_Bohr[3]; // x, y, and z-directional cell parameters, in atomic unit	

	// HOST_WORLD, GRID DATA
	int N[3];			 // Global 3-D grid dimensions
	MKL_Complex8 *cplx_gauss;	// (input) cellParams unit Bohr
	MKL_Complex8 *cplx_Ptr;

} compute_Models;

#endif

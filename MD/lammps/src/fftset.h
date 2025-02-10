#ifndef _FFTSET_H
#define _FFTSET_H

#include "mkl_types.h"
#include "tagmodels.h"
#include "computemodels.h"

int RealScaleImagZero(MKL_Complex8 *arr, int size, float scale);
int ComplexPointwiseMul(MKL_Complex8 *a, MKL_Complex8 *b, int size);

int gen_kernel(tag_Models **gk, double sigma);

int gaussR2K(compute_Models **compute_info);

int extfield_kernel(compute_Models **compute_info);

#endif

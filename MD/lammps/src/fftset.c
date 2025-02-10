#include "stdio.h"

#include "math.h"
#include "mkl.h"
#include "mkl_dfti.h"

#include "tagmodels.h"
#include "fftset.h"
/*
#ifdef _OPENMP
#include "omp.h"
#else
#define omp_get_thread_num() 0
#endif
*/
int RealScaleImagZero(MKL_Complex8 *arr, int size, float scale)
{
	/*
	int numThreads, threadID;

	#pragma omp parallel private(numThreads) private(threadID)
	{
		numThreads = omp_get_num_threads();
		threadID = omp_get_thread_num();

		for (int i = threadID; i < size; i += numThreads)
		{
			arr[i].real *= scale;
			arr[i].imag = 0.0;
		}
	}
	*/
	for (int i = 0; i < size; i++) {
		arr[i].real *= scale;
		arr[i].imag = 0.0;
	}
	
	return 0;
}

int ComplexPointwiseMul(MKL_Complex8 *a, MKL_Complex8 *b, int size)
{
	/*
	int numThreads, threadID;
	MKL_Complex8 c;

	#pragma omp parallel private(numThreads) private(threadID) private(c)
	{
		numThreads = omp_get_num_threads();
		threadID = omp_get_thread_num();

		for (int i = threadID; i < size; i += numThreads)
		{
			c.real = a[i].real * b[i].real - a[i].imag * b[i].imag;
			c.imag = a[i].real * b[i].imag + a[i].imag * b[i].real;
			b[i].real = c.real;
			b[i].imag = c.imag;
		}
	}
	*/
	MKL_Complex8 c;
	for (int i = 0; i < size; i++) {
		c.real = a[i].real * b[i].real - a[i].imag * b[i].imag;
		c.imag = a[i].real * b[i].imag + a[i].imag * b[i].real;
		b[i].real = c.real;
		b[i].imag = c.imag;		
	}
	
	
	return 0;
}

// gen_kernel
// input: Any of tag_Models' that can give units and size of gk. info.plrzblty[0] is a good choice. cellParams unit Bohr, pValue unit Ha/eBohr
// sigma: Gaussian kernel sigma, unit Bohr
// gk: Normalized Gaussian kernel, cellParams unit Bohr. gk MUST be already allocated with its memory, before this function call!
int gen_kernel(tag_Models **gk, double sigma)
{
	int i, j, k, rm_index;				 //rm_index is row-major order index.
	double s_sigma = sigma * sigma;		 // Prefix s_ stands for "squared"
	double r_exp_const = -0.5 / s_sigma; // Real space exponential constant
	double di, dj, dk, s_di, s_dj, s_dk; //di, dj, dk are displacements in x, y, and z directions.
	double i_test_integral;
	double grid_spacing[3], s_grid_spacing[3];
	double test_integral = 0.0;
	//double retest_integral = 0.0;

	for (i = 0; i < 3; i++)
	{
		grid_spacing[i] = (*gk)->cellParams[i] / (*gk)->pGrid[i];
		s_grid_spacing[i] = grid_spacing[i] * grid_spacing[i];
	}

	// Filling a real-space Gaussian distribution in gk.
	for (i = 0; i < (*gk)->pGrid[0]; i++)
	{
		for (j = 0; j < (*gk)->pGrid[1]; j++)
		{
			for (k = 0; k < (*gk)->pGrid[2]; k++)
			{
				rm_index = k + j * (*gk)->pGrid[2] + i * (*gk)->pGrid[2] * (*gk)->pGrid[1];

				di = i;
				dj = j;
				dk = k;

				di -= (double)rint(di / (*gk)->pGrid[0]) * (*gk)->pGrid[0];
				dj -= (double)rint(dj / (*gk)->pGrid[1]) * (*gk)->pGrid[1];
				dk -= (double)rint(dk / (*gk)->pGrid[2]) * (*gk)->pGrid[2];

				s_di = di * di;
				s_dj = dj * dj;
				s_dk = dk * dk;

				(*gk)->pValue[rm_index] = exp(r_exp_const * (s_grid_spacing[0] * s_di + s_grid_spacing[1] * s_dj + s_grid_spacing[2] * s_dk));
				test_integral += (*gk)->pValue[rm_index];
			}
		}
	}

	// Re-normalization of 3-D Gaussian
	i_test_integral = 1.0 / test_integral;
	for (i = 0; i < (*gk)->pGrid[0] * (*gk)->pGrid[1] * (*gk)->pGrid[2]; i++)
	{
		(*gk)->pValue[i] *= i_test_integral;
		//retest_integral += (*gk)->pValue[i];
	}

	return 0;
}

// gaussR2K
// Forward cuFFT for devcplx_gauss
int gaussR2K(compute_Models **compute_info)
{
	DFTI_DESCRIPTOR_HANDLE plan_gauss = NULL;
	MKL_LONG status;
	MKL_LONG dim_sizes[3] = {(*compute_info)->N[0], (*compute_info)->N[1], (*compute_info)->N[2]};

	/* Create a 3D FFT plan. */
	status = DftiCreateDescriptor(&plan_gauss, DFTI_SINGLE, DFTI_COMPLEX, 3, dim_sizes);
	status = DftiCommitDescriptor(plan_gauss);

	// Forward FFT
	status = DftiComputeForward(plan_gauss, (*compute_info)->cplx_gauss);

	// delete FFT plan
	status = DftiFreeDescriptor(&plan_gauss);

	return 0;
}

// extfield_kernel
// Gaussian operation in k-space
int extfield_kernel(compute_Models **compute_info)
{
	DFTI_DESCRIPTOR_HANDLE plan_extfield_kernel = NULL;
	MKL_LONG status;
	MKL_LONG dim_sizes[3] = {(*compute_info)->N[0], (*compute_info)->N[1], (*compute_info)->N[2]};
	
	int size = (*compute_info)->N[0] * (*compute_info)->N[1] * (*compute_info)->N[2];
	float norm_const = 1.0 / (float)(size);

	/* Create a 3D FFT plan. */
	status = DftiCreateDescriptor(&plan_extfield_kernel, DFTI_SINGLE, DFTI_COMPLEX, 3, dim_sizes);
	status = DftiCommitDescriptor(plan_extfield_kernel);

	// Forward FFT
	status = DftiComputeForward(plan_extfield_kernel, (*compute_info)->cplx_Ptr);

	// Gaussian operation in k-space
	ComplexPointwiseMul((*compute_info)->cplx_gauss, (*compute_info)->cplx_Ptr, \
	(*compute_info)->N[0] * (*compute_info)->N[1] * (*compute_info)->N[2]);

	// Backward FFT
	status = DftiComputeBackward(plan_extfield_kernel, (*compute_info)->cplx_Ptr);

	// delete FFT plan
	status = DftiFreeDescriptor(&plan_extfield_kernel);

	// Normalizing FFT arrays
	RealScaleImagZero((*compute_info)->cplx_Ptr, (*compute_info)->N[0] * (*compute_info)->N[1] * (*compute_info)->N[2], norm_const);

	return 0;
}

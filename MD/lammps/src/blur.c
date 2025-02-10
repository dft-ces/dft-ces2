#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "mkl_types.h"
#include "unistd.h"

#include "tagmodels.h"
#include "computemodels.h"
#include "fftset.h"
#include "blur.h"

/* ---------------------------------------------------------------------- */

int main(int argc, char **argv)
{
	int i, j, k, dir, rm_index;
	char potFile[200];
	double sigma;
	struct timespec stampParsing1, stampParsing2, \
		stampCopybyvalue1, stampCopybyvalue2, \
		stampFFT1, stampFFT2, \
		stampSave1, stampSave2;
	double durationParsing, durationCopybyvalue, durationFFT, durationSave;

	// DEBUG PART START
	/*
	double test_integral;
	double idV;
	*/
	// DEBUG PART END

	// Argument parsing
	if (argc == 3) {
		sprintf(potFile,"%s",argv[1]);
		if(access(potFile,F_OK) == -1) {
			printf("File does not exist.\n");
			return 1;
		}
		sigma = atof(argv[2]);
	}
	else {
		printf("Usage: ./blur.o cubefile sigma\n");
		return 1;
	}

	// tag_Models: a struct for cubefile data
	tag_Models *base, *gk;

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampParsing1);
	potParser(&base, potFile);
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampParsing2);
	durationParsing = (1000.0 * (double)stampParsing2.tv_sec + 1.0e-6 * stampParsing2.tv_nsec) - \
	(1000.0 * (double)stampParsing1.tv_sec + 1.0e-6 * stampParsing1.tv_nsec);

	potCopywoValue(&base, &gk);

	// Generating normalized real-space Gaussian centered at origin
	gen_kernel(&gk, sigma);
	printf("Gaussian sigma: %13.5E\n", sigma);

	// DEBUG PART START
	//saveCube(&gk, "kernel.cube");
	/*
	idV = (base->pGrid[0] * base->pGrid[1] * base->pGrid[2]) / (base->cellParams[0] * base->cellParams[1] * base->cellParams[2]);
	test_integral = 0.0;
	for (i = 0; i < base->pGrid[0] * base->pGrid[1] * base->pGrid[2]; i++) {
		test_integral += base->pValue[i];
	}
	test_integral *= idV;
	printf("3D Integral of the input grid: %13.5E\n", test_integral);
	*/
	// DEBUG PART END

	// compute_Models: a struct for cuFFT data
	compute_Models *compute_info;
	compute_info = (compute_Models *)malloc(sizeof(compute_Models));

	for (dir = 0; dir < 3; dir++) {
		compute_info->cell_Bohr[dir] = base->cellParams[dir];
		compute_info->N[dir] = base->pGrid[dir];
	}

	compute_info->cplx_gauss = (MKL_Complex8 *)malloc(sizeof(MKL_Complex8) * compute_info->N[0] * compute_info->N[1] * compute_info->N[2]);
	compute_info->cplx_Ptr = (MKL_Complex8 *)malloc(sizeof(MKL_Complex8) * compute_info->N[0] * compute_info->N[1] * compute_info->N[2]);

	// copying data into host pointers
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampCopybyvalue1);
	for (i = 0; i < compute_info->N[0]; i++) {
		for (j = 0; j < compute_info->N[1]; j++) {
			for (k = 0; k < compute_info->N[2]; k++) {
				rm_index = k + j * compute_info->N[2] + i * compute_info->N[2] * compute_info->N[1];
				compute_info->cplx_gauss[rm_index].real = (gk->pValue[rm_index]);
				compute_info->cplx_gauss[rm_index].imag = 0.0;				
				compute_info->cplx_Ptr[rm_index].real = base->pValue[rm_index];
				compute_info->cplx_Ptr[rm_index].imag = 0.0;
			}
		}
	}
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampCopybyvalue2);
	durationCopybyvalue = (1000.0 * (double)stampCopybyvalue2.tv_sec + 1.0e-6 * stampCopybyvalue2.tv_nsec) - \
	(1000.0 * (double)stampCopybyvalue1.tv_sec + 1.0e-6 * stampCopybyvalue1.tv_nsec);

	// gaussR2K
	// Forward discrete Fourier transform for devcplx_gauss
	gaussR2K(&compute_info);

	// extfield_kernel
	// Gaussian operation in k-space
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampFFT1);
	extfield_kernel(&compute_info);
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampFFT2);
	durationFFT = (1000.0 * (double)stampFFT2.tv_sec + 1.0e-6 * stampFFT2.tv_nsec) - \
	(1000.0 * (double)stampFFT1.tv_sec + 1.0e-6 * stampFFT1.tv_nsec);

	// copying data from the host pointer to the tag_Models type pointer
	for (i = 0; i < compute_info->N[0]; i++) {
		for (j = 0; j < compute_info->N[1]; j++) {
			for (k = 0; k < compute_info->N[2]; k++) {
				rm_index = k + j * compute_info->N[2] + i * compute_info->N[2] * compute_info->N[1];
				base->pValue[rm_index] = compute_info->cplx_Ptr[rm_index].real;
			}
		}
	}

	// DEBUG PART START
	/*
	test_integral = 0.0;
	for (i = 0; i < base->pGrid[0] * base->pGrid[1] * base->pGrid[2]; i++) {
		test_integral += base->pValue[i];
	}
	test_integral *= idV;
	printf("3D Integral of the blurred grid: %13.5E\n", test_integral);
	*/
	// DEBUG PART END

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampSave1);
	saveCube(&base, "blurred.cube");
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &stampSave2);
	durationSave = (1000.0 * (double)stampSave2.tv_sec + 1.0e-6 * stampSave2.tv_nsec) - \
	(1000.0 * (double)stampSave1.tv_sec + 1.0e-6 * stampSave1.tv_nsec);

	printf("\n");
	printf("Time taken to parse a cube file: %f ms\n", durationParsing);
	printf("Time taken to copy grid data into a complex array: %f ms\n", durationCopybyvalue);
	printf("Time taken to perform the gaussian kernel operation (2 FFTs): %f ms\n", durationFFT);
	printf("Time taken to save a cube file: %f ms\n", durationSave);

	free(compute_info->cplx_gauss);
	free(compute_info->cplx_Ptr);
	free(compute_info);
	potFree(&base);
	potFree(&gk);

	return 0;
}

// potParser
// target cube file: potential, cellParams unit Bohr, pValue in Hartree atomic units,
// returned tag_Models: potential, cellParams unit Bohr, pValue in Hartree atomic units
int potParser(tag_Models **output, char *fname)
{ // LOCPOT parser
	FILE *fp;
	char str[2048], *str_ptr;
	int nline = 0, cnt = 0;

	*output = (tag_Models *)malloc(sizeof(tag_Models));

	printf("###  Grid file info : %s\n", fname);

	fp = fopen(fname, "r");
	(*output)->nAtoms = 0;
	(*output)->cellParams[0] = (*output)->cellParams[1] = (*output)->cellParams[2] = (double)0;

	while (fgets(str, sizeof(str), fp) != NULL)
	{
		int i, j, *temp, k;
		double a, b, c, d;
		nline++;

		if (nline == 3)
		{
			sscanf(str, "%d %*f %*f %*f", &(*output)->nAtoms);
			(*output)->atomN = (int *)malloc((*output)->nAtoms * sizeof(int));
			(*output)->posX = (double *)malloc((*output)->nAtoms * sizeof(double));
			(*output)->posY = (double *)malloc((*output)->nAtoms * sizeof(double));
			(*output)->posZ = (double *)malloc((*output)->nAtoms * sizeof(double));
		}

		if (nline == 4)
		{
			sscanf(str, "%d %lf %*f %*f", &(*output)->pGrid[0], &(*output)->cellParams[0]);
			(*output)->cellParams[0] = (*output)->pGrid[0] * (*output)->cellParams[0];
		}

		if (nline == 5)
		{
			sscanf(str, "%d %*f %lf %*f", &(*output)->pGrid[1], &(*output)->cellParams[1]);
			(*output)->cellParams[1] = (*output)->pGrid[1] * (*output)->cellParams[1];
		}

		if (nline == 6)
		{
			sscanf(str, "%d %*f %*f %lf", &(*output)->pGrid[2], &(*output)->cellParams[2]);
			(*output)->cellParams[2] = (*output)->pGrid[2] * (*output)->cellParams[2];
		}

		if (nline >= 7 && nline < 7 + (*output)->nAtoms)
		{
			sscanf(str, "%d %*f %lf %lf %lf", &(*output)->atomN[nline - 7], &(*output)->posX[nline - 7], &(*output)->posY[nline - 7], &(*output)->posZ[nline - 7]);
		}

		if (nline == (*output)->nAtoms + 6)
		{
			printf("### Cell parameters(A): %f %f %f\n", (*output)->cellParams[0], (*output)->cellParams[1], (*output)->cellParams[2]);
			printf("### Number of Atoms: %d\n", (*output)->nAtoms);
			printf("### Grid Dimension: %d %d %d\n", (*output)->pGrid[0], (*output)->pGrid[1], (*output)->pGrid[2]);
			(*output)->pValue = (double *)malloc((*output)->pGrid[0] * (*output)->pGrid[1] * (*output)->pGrid[2] * sizeof(double));
		}

		if (nline >= (*output)->nAtoms + 7 && cnt < (*output)->pGrid[0] * (*output)->pGrid[1] * (*output)->pGrid[2])
		{
			str_ptr = strtok(str, " ");
			for (; str_ptr != NULL; cnt++)
			{
				sscanf(str_ptr, "%lf", &(*output)->pValue[cnt]);
				(*output)->pValue[cnt] = (*output)->pValue[cnt];
				str_ptr = strtok(NULL, " ");
			}
		}

	} //while end

	return 0;
}

// saveCube: saves a tag_Models input into a cube file
int saveCube(tag_Models **input, char *fname)
{
	FILE *fp;
	int i, j, k, cnt;
	fp = fopen(fname, "w");

	fprintf(fp, " Gaussian Blurred Cubefile\n");
	fprintf(fp, " lattice unit Bohr\n");
	fprintf(fp, "% 5d    0.000000    0.000000    0.000000\n", (*input)->nAtoms);
	fprintf(fp, "% 5d% 12.6lf    0.000000    0.000000\n", (*input)->pGrid[0], (*input)->cellParams[0] / (*input)->pGrid[0]);
	fprintf(fp, "% 5d    0.000000% 12.6lf    0.000000\n", (*input)->pGrid[1], (*input)->cellParams[1] / (*input)->pGrid[1]);
	fprintf(fp, "% 5d    0.000000    0.000000% 12.6lf\n", (*input)->pGrid[2], (*input)->cellParams[2] / (*input)->pGrid[2]);
	for (i = 0; i < (*input)->nAtoms; i++)
	{
		fprintf(fp, "% 5d% 12.6lf% 12.6lf% 12.6lf% 12.6lf\n", (*input)->atomN[i], (double)(*input)->atomN[i], (*input)->posX[i], (*input)->posY[i], (*input)->posZ[i]);
	}
	cnt = 0;
	for (i = 0; i < (*input)->pGrid[0]; i++)
	{
		for (j = 0; j < (*input)->pGrid[1]; j++)
		{
			for (k = 0; k < (*input)->pGrid[2]; k++)
			{
				fprintf(fp, "% 13.5lE", ((*input)->pValue[k + j * (*input)->pGrid[2] + i * (*input)->pGrid[2] * (*input)->pGrid[1]]));
				if (cnt % 6 == 5)
					fprintf(fp, "\n");
				cnt++;
			}
		}
	}
	fclose(fp);
	printf("### New cubefile has been saved as %s\n", fname);

	return 0;
}

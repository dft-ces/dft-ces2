#include "stdlib.h"
#include "tagmodels.h"

int potCopywoValue(tag_Models **in, tag_Models **out)
{
	int i;
	*out = (tag_Models *)malloc(sizeof(tag_Models));
	(*out)->nAtoms = (*in)->nAtoms;
	(*out)->cellParams[0] = (*in)->cellParams[0];
	(*out)->cellParams[1] = (*in)->cellParams[1];
	(*out)->cellParams[2] = (*in)->cellParams[2];
	(*out)->pGrid[0] = (*in)->pGrid[0];
	(*out)->pGrid[1] = (*in)->pGrid[1];
	(*out)->pGrid[2] = (*in)->pGrid[2];
	(*out)->atomN = (int *)malloc(sizeof(int) * (*in)->nAtoms);
	(*out)->posX = (double *)malloc(sizeof(double) * (*in)->nAtoms);
	(*out)->posY = (double *)malloc(sizeof(double) * (*in)->nAtoms);
	(*out)->posZ = (double *)malloc(sizeof(double) * (*in)->nAtoms);

	for (i = 0; i < (*in)->nAtoms; i++)
	{
		(*out)->atomN[i] = (*in)->atomN[i];
		(*out)->posX[i] = (*in)->posX[i];
		(*out)->posY[i] = (*in)->posY[i];
		(*out)->posZ[i] = (*in)->posZ[i];
	}

	(*out)->pValue = (double *)malloc((*in)->pGrid[0] * (*in)->pGrid[1] * (*in)->pGrid[2] * sizeof(double));

	return 0;
}

int potFree(tag_Models **in)
{
	free((*in)->pValue);
	free((*in)->atomN);
	free((*in)->posX);
	free((*in)->posY);
	free((*in)->posZ);
	free(*in);
}

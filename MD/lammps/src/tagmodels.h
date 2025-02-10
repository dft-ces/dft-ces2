#ifndef _TAGMODELS_H
#define _TAGMODELS_H

typedef struct
{ //LOCPOT parsing data structure
	double cellParams[3];
	int nAtoms;
	int pGrid[3]; //grid x, y, z
	int *atomN;
	double *posX;
	double *posY;
	double *posZ;
	double *pValue;
} tag_Models;

int potCopywoValue(tag_Models **in, tag_Models **out);
int potFree(tag_Models **in);

#endif

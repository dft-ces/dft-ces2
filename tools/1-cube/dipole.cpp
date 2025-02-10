#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string.h>
#include <unistd.h>

typedef struct {
  double cellParams[3];
  int nAtoms;
  int pGrid[3];
  int *atomN;
  double *posX;
  double *posY;
  double *posZ;
  double *pValue;
} tag_Models;

tag_Models chgParser(char *fname);
void saveCube(tag_Models input, char *fname);
double dipx,dipy,dipz;
main(int argc, char *argv[]) {
  int i, j, k, dir, position;
  char chgFile[200];
  double dipole, ratio;
  // arguments parsing
  if (argc == 2) {
    sprintf(chgFile, "%s", argv[1]);
    if (access(chgFile, F_OK) == -1) {
      printf("File does not exist.\n");
      return 0;
    }
  } else {
    printf("Usage: ./dipc chgden.cube \n");
    return 0;
  }

  tag_Models baseModel;
  baseModel = chgParser(chgFile);
  int cnt = 3;
  double gVol = baseModel.cellParams[0] * baseModel.cellParams[1] * baseModel.cellParams[2] / (baseModel.pGrid[0] * baseModel.pGrid[1] * baseModel.pGrid[2]);

  if (cnt == 3) {
    dipole = 0;
    position = int(baseModel.pGrid[2]);
    for (i = 0; i < baseModel.pGrid[0]; i++) {
      for (j = 0; j < baseModel.pGrid[1]; j++) {
        for (k = 0; k < baseModel.pGrid[2]; k++) {
          if (k < position) {
            dipole += baseModel.pValue[k + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] * k * baseModel.cellParams[2] / baseModel.pGrid[2];
          } else {
            dipole += baseModel.pValue[k + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] * (k - baseModel.pGrid[2]) * baseModel.cellParams[2] / baseModel.pGrid[2];
          }
        }
      }
    }
    //    for (i = 0; i < baseModel.pGrid[0]; i++) {
    //      for (j = 0; j < baseModel.pGrid[1]; j++) {
    //        baseModel.pValue[position + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] += dipole / baseModel.pGrid[0] / baseModel.pGrid[1];
    //        baseModel.pValue[(position + 1) + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] -= dipole / baseModel.pGrid[0] / baseModel.pGrid[1];
    //      }
    //    }
		dipz = dipole * gVol / 0.3934303;
    printf("dipole moment in z direction is %f D \n", dipz);
  }
  cnt = 2;
  if (cnt == 2) {
    dipole = 0;
    position = int(baseModel.pGrid[1]);
    for (i = 0; i < baseModel.pGrid[0]; i++) {
      for (k = 0; j < baseModel.pGrid[2]; k++) {
        for (j = 0; j < baseModel.pGrid[1]; j++) {
          if (j < position) {
            dipole += baseModel.pValue[k + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] * j * baseModel.cellParams[1] / baseModel.pGrid[1];
          } else {
            dipole += baseModel.pValue[k + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] * (j - baseModel.pGrid[1]) * baseModel.cellParams[1] / baseModel.pGrid[1];
          }
        }
      }
    }
    // for (i = 0; i < baseModel.pGrid[0]; i++) {
    //   for (k = 0; k < baseModel.pGrid[2]; k++) {
    //     baseModel.pValue[k + position * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] += dipole / baseModel.pGrid[0] / baseModel.pGrid[2];
    //     baseModel.pValue[k + (position + 1) * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] -= dipole / baseModel.pGrid[0] / baseModel.pGrid[2];
    //   }
    // }
		dipy = dipole * gVol / 0.3934303;
    printf("dipole moment in y direction is %f D \n", dipy);
  }
  cnt = 1;
  if (cnt == 1) {
    dipole = 0;
    position = int(baseModel.pGrid[0]);
    for (j = 0; j < baseModel.pGrid[1]; j++) {
      for (k = 0; k < baseModel.pGrid[2]; k++) {
        for (i = 0; i < baseModel.pGrid[0]; i++) {
          if (i < position) {
            dipole += baseModel.pValue[k + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] * i * baseModel.cellParams[0] / baseModel.pGrid[0];
          } else {
            dipole += baseModel.pValue[k + j * baseModel.pGrid[2] + i * baseModel.pGrid[2] * baseModel.pGrid[1]] * (i - baseModel.pGrid[0]) * baseModel.cellParams[0] / baseModel.pGrid[0];
          }
        }
      }
    }
    // for (j = 0; j < baseModel.pGrid[1]; j++) {
    //   for (k = 0; k < baseModel.pGrid[2]; k++) {
    //     baseModel.pValue[k + j * baseModel.pGrid[2] + position * baseModel.pGrid[2] * baseModel.pGrid[1]] += dipole / baseModel.pGrid[2] / baseModel.pGrid[1];
    //     baseModel.pValue[k + j * baseModel.pGrid[2] + (position + 1) * baseModel.pGrid[2] * baseModel.pGrid[1]] -= dipole / baseModel.pGrid[2] / baseModel.pGrid[1];
    //   }
    // }
		dipx = dipole * gVol / 0.3934303;
    printf("dipole moment in x direction is %f D \n", dipx);
    printf("Magnitude of the dipole moment is %f D \n", sqrt(dipx*dipx+dipy*dipy+dipz*dipz));
  }
  //  saveCube(baseModel, "dipc.cube");
} // main END

tag_Models chgParser(char *fname) { // charge density file parsing using atomic units
  FILE *fp;
  char str[2048], *str_ptr;
  int flag = 0, nline = 0, cnt = 0, nAtoms = 0;
  tag_Models output;

  fp = fopen(fname, "r");
  output.nAtoms = 0;
  output.cellParams[0] = output.cellParams[1] = output.cellParams[2] = (double)0;

  while (fgets(str, sizeof(str), fp) != NULL) {
    int i, j, *temp;
    double a, b, c;
    nline++;

    if (nline == 3) {
      sscanf(str, "%d %*f %*f %*f", &output.nAtoms);
      output.atomN = (int *)malloc(output.nAtoms * sizeof(int));
      output.posX = (double *)malloc(output.nAtoms * sizeof(double));
      output.posY = (double *)malloc(output.nAtoms * sizeof(double));
      output.posZ = (double *)malloc(output.nAtoms * sizeof(double));
    }

    if (nline == 4) {
      sscanf(str, "%d %lf %*f %*f", &output.pGrid[0], &output.cellParams[0]);
      output.cellParams[0] = output.pGrid[0] * output.cellParams[0];
    }

    if (nline == 5) {
      sscanf(str, "%d %*f %lf %*f", &output.pGrid[1], &output.cellParams[1]);
      output.cellParams[1] = output.pGrid[1] * output.cellParams[1];
    }

    if (nline == 6) {
      sscanf(str, "%d %*f %*f %lf", &output.pGrid[2], &output.cellParams[2]);
      output.cellParams[2] = output.pGrid[2] * output.cellParams[2];
    }

    if (nline >= 7 && nline < 7 + output.nAtoms) {
      sscanf(str, "%d %*f %lf %lf %lf", &output.atomN[nline - 7], &output.posX[nline - 7], &output.posY[nline - 7], &output.posZ[nline - 7]);
    }

    if (nline == output.nAtoms + 6) {
      printf("### Cell parameters(Bohr): %f %f %f\n", output.cellParams[0], output.cellParams[1], output.cellParams[2]);
      printf("### Number of Atoms: %d\n", output.nAtoms);
      printf("### Grid Dimension: %d %d %d\n", output.pGrid[0], output.pGrid[1], output.pGrid[2]);
      output.pValue = (double *)malloc(output.pGrid[0] * output.pGrid[1] * output.pGrid[2] * sizeof(double));
    }

    if (nline >= output.nAtoms + 7 && cnt < output.pGrid[0] * output.pGrid[1] * output.pGrid[2]) {
      str_ptr = strtok(str, " ");
      for (; str_ptr != NULL; cnt++) {
        sscanf(str_ptr, "%lf", &output.pValue[cnt]);
        str_ptr = strtok(NULL, " ");
      }
    }

  } // while end
  return output;
}

void saveCube(tag_Models input, char *fname) {
  FILE *fp;
  int i, j, k, cnt;
  fp = fopen(fname, "w");

  fprintf(fp, " Cubefile created from chg2pot\n");
  fprintf(fp, " ES potential, lattice unit Bohr, grid unit Ry\n");
  fprintf(fp, "% 5d    0.000000    0.000000    0.000000\n", input.nAtoms);
  fprintf(fp, "% 5d% 12.6lf    0.000000    0.000000\n", input.pGrid[0], input.cellParams[0] / input.pGrid[0]);
  fprintf(fp, "% 5d    0.000000% 12.6lf    0.000000\n", input.pGrid[1], input.cellParams[1] / input.pGrid[1]);
  fprintf(fp, "% 5d    0.000000    0.000000% 12.6lf\n", input.pGrid[2], input.cellParams[2] / input.pGrid[2]);
  for (i = 0; i < input.nAtoms; i++) {
    fprintf(fp, "% 5d% 12.6lf% 12.6lf% 12.6lf% 12.6lf\n", input.atomN[i], (double)input.atomN[i], input.posX[i], input.posY[i], input.posZ[i]);
  }
  cnt = 0;
  for (i = 0; i < input.pGrid[0]; i++) {
    for (j = 0; j < input.pGrid[1]; j++) {
      for (k = 0; k < input.pGrid[2]; k++) {
        fprintf(fp, "% 13.5lE", (input.pValue[k + j * input.pGrid[2] + i * input.pGrid[2] * input.pGrid[1]]));
        if (cnt % 6 == 5)
          fprintf(fp, "\n");
        cnt++;
      }
    }
  }
  fclose(fp);
  printf("### Potential was saved to %s\n", fname);
}

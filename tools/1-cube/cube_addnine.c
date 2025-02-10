#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

tag_Models potParser(char *fname);
void saveCube(tag_Models input, char *fname);

main(int argc, char *argv[]) {
  int i, j, k, dir;
  char potFile1[200];
  char potFile2[200];
  char potFile3[200];
  char potFile4[200];
  char potFile5[200];
  char potFile6[200];
  char potFile7[200];
  char potFile8[200];
  char potFile9[200];
  // char potFile10[200];
  // char potFile11[200];
  //  arguments parsing
  if (argc == 10) {
    sprintf(potFile1, "%s", argv[1]);
    if (access(potFile1, F_OK) == -1) {
      printf("File does not exist.\n");
      return 0;
    }
    sprintf(potFile2, "%s", argv[2]);
    if (access(potFile2, F_OK) == -1) {
      printf("File does not exist.\n");
      return 0;
    }
    sprintf(potFile3, "%s", argv[3]);
    if (access(potFile2, F_OK) == -1) {
      printf("File does not exist.\n");
      return 0;
    }
    sprintf(potFile4, "%s", argv[4]);
    sprintf(potFile5, "%s", argv[5]);
    sprintf(potFile6, "%s", argv[6]);
    sprintf(potFile7, "%s", argv[7]);
    sprintf(potFile8, "%s", argv[8]);
    sprintf(potFile9, "%s", argv[9]);
    // sprintf(potFile10, "%s", argv[10]);
    // sprintf(potFile11, "%s", argv[11]);
    if (access(potFile2, F_OK) == -1) {
      printf("File does not exist.\n");
      return 0;
    }
  } else {
    printf("for cube1+cube2+cube3+cube4+cube..n, Usage: ./cube_add cube1 cube2 cube3 cube4 .. cuben\n");
    return 0;
  }

  tag_Models base1, base2, base3, base4, base5, base6, base7, base8, base9;
  base1 = potParser(potFile1);
  base2 = potParser(potFile2);
  base3 = potParser(potFile3);
  base4 = potParser(potFile4);
  base5 = potParser(potFile5);
  base6 = potParser(potFile6);
  base7 = potParser(potFile7);
  base8 = potParser(potFile8);
  base9 = potParser(potFile9);
  // base10 = potParser(potFile10);
  // base11 = potParser(potFile11);

  for (i = 0; i < base1.pGrid[0]; i++) {
    for (j = 0; j < base1.pGrid[1]; j++) {
      for (k = 0; k < base1.pGrid[2]; k++) {
        base1.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] += base2.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base3.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base4.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base5.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base6.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base7.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base8.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]] + base9.pValue[k + j * base1.pGrid[2] + i * base1.pGrid[2] * base1.pGrid[1]];
      }
    }
  }

  saveCube(base1, "add.cube");
} // main END

tag_Models potParser(char *fname) { // LOCPOT parser
  FILE *fp;
  char str[2048], *str_ptr;
  int flag = 0, nline = 0, cnt = 0, nAtoms = 0;
  tag_Models output;

  fp = fopen(fname, "r");
  output.nAtoms = 0;
  output.cellParams[0] = output.cellParams[1] = output.cellParams[2] = (double)0;

  while (fgets(str, sizeof(str), fp) != NULL) {
    int i, j, *temp, k;
    double a, b, c, d;
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

  fprintf(fp, " Cubefile rescaled\n");
  fprintf(fp, " lattice unit Bohr\n");
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
  printf("### New cubefile has been saved as %s\n", fname);
}

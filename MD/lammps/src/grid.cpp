/* ----------------------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#include "grid.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_gridforce.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "blur.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

Grid::Grid(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp) {
  int me = comm->me;
  int i;

  natoms = 0;
  atomns = NULL;
  basis = NULL;
  gvin = NULL;
  gvout = NULL;
  gvout_all = NULL;
  cvin = NULL;
  savedcube = NULL;
  inputrhocube = NULL;

  if (narg < 4)
    error->all(FLERR, "Illegal grid command"); // e.g. n+1(# of input cubes) n'+5(# of output cubes) pot.cube|None rho_H.cube|None rho_O.cube hy.cube ox.cube repA.cube dipx.cube dipy.cube dipz.cube

  if (strcmp(arg[2], "None") != 0) {
    // read header until atom information
    read_header(arg[2]);

    if (gx[1] + gx[2] + gy[0] + gy[2] + gz[0] + gz[1] != 0)
      error->all(FLERR, "Grid file is not orthorhombic"); // gx[0], gy[1], gz[2] are not zero.
    if (me == 0) {
      printf("Grid file will be parsed: %s\n", arg[2]);
      printf("DFT-CES: # of atoms in grid: %d\n", natoms);
      printf("DFT-CES: # of grid points: %d %d %d\n", gnx, gny, gnz);
      printf("DFT-CES: grid spacing: %f %f %f Angs\n", gx[0], gy[1], gz[2]);
    }
    memory->grow(atomns, natoms, "grid:atomns");
    memory->grow(basis, natoms, 3, "grid:basis");
    memory->grow(gvin, gnx * gny * gnz, "grid:gvin");
    // read atom and grid value information
    read_content(arg[2]);
  } else {
    return;
  }

  nrhocubes = atoi(arg[0]) - 1;                    // arg[0] pot.cube rho_x.cube : n(rho cubes)+1(pot.cube)
                                                   // parsing cube file of charge density.
  if (strcmp(arg[3], "None") != 0) {               // int int pot.cube "the first rho cube"  at arg[3].
                                                   // read header until atom information
    inputrhocube = new char *[nrhocubes];          // cubefile for input QM rho cubes
    for (int temp = 0; temp < nrhocubes; temp++) { // usually "nrhocubes" is 2.
      inputrhocube[temp] = (char *)"empty";
    }
    for (int temp = 0; temp < nrhocubes; temp++) {
      int n = strlen(arg[temp + 3] + 1);
      inputrhocube[temp] = new char[n];
      strcpy(inputrhocube[temp], arg[temp + 3]);
      if (me == 0)
        printf("Grid file will be parsed: %s\n", inputrhocube[temp]);
    }
    memory->grow(cvin, gnx * gny * gnz * nrhocubes, "grid:cvin");

    // read atom and grid value information
    read_chd(inputrhocube);
  } else {
    if (me == 0)
      printf("charge density Grid file will not be parsed. \n");
    return;
  }
  savetag = 1; // default is save
  if (strcmp(arg[1], "0") == 0) {
    savetag = 0; // don't save
  }
  if (savetag == 0) {
    if (me == 0)
      printf("Grid file will not be saved\n");
  } else { // Grid file will be saved.

    ncubes = atoi(arg[1]);
    memory->grow(gvout, gnx * gny * gnz * (atoi(arg[1])), "grid:gvout"); // atoi(arg[0]) - 1 + 3 + tip4p_CES, +3 for dipx, dipy, dipz
    if (me == 0)
      memory->grow(gvout_all, gnx * gny * gnz * (atoi(arg[1])), "grid:gvout_all");

    savedcube = new char *[ncubes];             // cubefile for save
    for (int temp = 0; temp < ncubes; temp++) { // usually "ncubes" is 6.
      savedcube[temp] = (char *)"empty";
    }
    for (int temp = 0; temp < ncubes; temp++) {
      int n = strlen(arg[temp + nrhocubes + 3] + 1); // int(nrhocubes+1) int nrhocubes+1 "the first save cube" at arg[nrhocubes+3].
      savedcube[temp] = new char[n];
      strcpy(savedcube[temp], arg[temp + nrhocubes + 3]);
      if (me == 0)
        printf("grid will be saved in %s\n", savedcube[temp]);
    }
  }
}
/* ---------------------------------------------------------------------- */

Grid::~Grid() {
  memory->destroy(atomns);
  memory->destroy(basis);
  memory->destroy(gvin);
  memory->destroy(gvout);
  memory->destroy(gvout_all);
  memory->destroy(cvin);
}

/* ---------------------------------------------------------------------- */

void Grid::read_header(char *filename) {
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];
  double temp[4];

  if (me == 0) {
    fptr = fopen(filename, "r");
    if (fptr == NULL) {
      char str[128];
      sprintf(str, "Cannot open grid file %s", filename);
      error->one(FLERR, str);
    }

    fgets(line, MAXLINE, fptr);
    fgets(line, MAXLINE, fptr); // skip first two rows
    fgets(line, MAXLINE, fptr);
    sscanf(line, "%d %lf %lf %lf", &natoms, &temp[0], &temp[1], &temp[2]);
    fgets(line, MAXLINE, fptr);
    sscanf(line, "%d %lf %lf %lf", &gnx, &gx[0], &gx[1], &gx[2]);
    gx[0] *= 0.52917721; // bohr to angs
    fgets(line, MAXLINE, fptr);
    sscanf(line, "%d %lf %lf %lf", &gny, &gy[0], &gy[1], &gy[2]);
    gy[1] *= 0.52917721;
    fgets(line, MAXLINE, fptr);
    sscanf(line, "%d %lf %lf %lf", &gnz, &gz[0], &gz[1], &gz[2]);
    gz[2] *= 0.52917721;

    fclose(fptr);
  }

  MPI_Bcast(&natoms, 1, MPI_INT, 0, world);
  MPI_Bcast(&gnx, 1, MPI_INT, 0, world);
  MPI_Bcast(&gny, 1, MPI_INT, 0, world);
  MPI_Bcast(&gnz, 1, MPI_INT, 0, world);
  MPI_Bcast(&gx[0], 3, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gy[0], 3, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gz[0], 3, MPI_DOUBLE, 0, world);
}

void Grid::read_content(char *filename) {
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE], *str_ptr;
  double temp;
  int i, cnt;

  if (me == 0) {
    fptr = fopen(filename, "r");
    for (i = 0; i < 6; i++)
      fgets(line, MAXLINE, fptr);
    for (i = 0; i < natoms; i++) {
      fgets(line, MAXLINE, fptr);
      sscanf(line, "%d %lf %lf %lf %lf", &atomns[i], &temp, &basis[i][0], &basis[i][1], &basis[i][2]);
      basis[i][0] *= 0.52917721;
      basis[i][1] *= 0.52917721;
      basis[i][2] *= 0.52917721;
    }
    cnt = 0;
    while (fgets(line, MAXLINE, fptr) != NULL) {
      str_ptr = strtok(line, " ");
      for (; str_ptr != NULL; cnt++) {
        sscanf(str_ptr, "%lf", &gvin[cnt]);
        gvin[cnt] *= 13.60569253; // Ry to eV
        str_ptr = strtok(NULL, " ");
      }
    }
    fclose(fptr);
    printf("\n#######\nDFT-CES: %s has been parsed\n\n", filename);
  }

  MPI_Bcast(&atomns[0], natoms, MPI_INT, 0, world);
  MPI_Bcast(&basis[0][0], 3 * natoms, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gvin[0], gnx * gny * gnz, MPI_DOUBLE, 0, world);
}

void Grid::read_chd(char **filename) // jth
{
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE], *str_ptr;
  int i, cnt;
  int block = gnx * gny * gnz;
  if (me == 0) {
    for (int temp = 0; temp < nrhocubes; temp++) { // g'daymate
      fptr = fopen(filename[temp], "r");
      if (fptr == NULL) {
        char str[128];
        sprintf(str, "Revised. Cannot open grid file %s", filename[temp]);
        error->one(FLERR, str);
      }
      for (i = 0; i < 6; i++)
        fgets(line, MAXLINE, fptr);
      for (i = 0; i < natoms; i++) {
        fgets(line, MAXLINE, fptr);
      }
      cnt = 0;
      while (fgets(line, MAXLINE, fptr) != NULL) {
        str_ptr = strtok(line, " ");
        for (; str_ptr != NULL; cnt++) {
          sscanf(str_ptr, "%lf", &cvin[block * temp + cnt]);
          //		else if(temp ==1) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          //		else if(temp ==2) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          //		else if(temp ==3) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          //		else if(temp ==4) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          //		else if(temp ==5) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          //		else if(temp ==6) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          //		else if(temp ==7) sscanf(str_ptr,"%lf",&cvin[temp*gnx*gny*gnz+cnt]);
          str_ptr = strtok(NULL, " ");
        }
      }
      fclose(fptr);
      printf("\n#######\nDFT-CES: %s has been parsed\n\n", filename[temp]);
    }
  } // open me == 0
  MPI_Bcast(&cvin[0], nrhocubes * gnx * gny * gnz, MPI_DOUBLE, 0, world);
  //  for (int temp=0;temp<nrhocubes;temp++){
  //	  if(temp ==0) MPI_Bcast(&cvin[0],gnx*gny*gnz,MPI_DOUBLE,0,world);
  //	  else if(temp ==1)  MPI_Bcast(&cvin1[0],gnx*gny*gnz,MPI_DOUBLE,0,world);
  //          else if(temp ==2)  MPI_Bcast(&cvin2[0],gnx*gny*gnz,MPI_DOUBLE,0,world);
  //          else if(temp ==3)  MPI_Bcast(&cvin3[0],gnx*gny*gnz,MPI_DOUBLE,0,world);
  //          else if(temp ==4)  MPI_Bcast(&cvin4[0],gnx*gny*gnz,MPI_DOUBLE,0,world);
  //  }
}

void Grid::save_grid(char **filename, int nsteps) {
  int me = comm->me;
  FILE *fptr;
  int i, j, k, cnt;
  int block;
  block = gnx * gny * gnz;
  MPI_Reduce(gvout, gvout_all, block * ncubes, MPI_DOUBLE, MPI_SUM, 0, world);
  for (int temp = 0; temp < ncubes; temp++) {
    //    if(temp ==0) MPI_Reduce(gvout0,gvout_all0,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);
    //    else if(temp ==1) MPI_Reduce(gvout1,gvout_all1,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);
    //    else if(temp ==2) MPI_Reduce(gvout2,gvout_all2,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);
    //    else if(temp ==3) MPI_Reduce(gvout3,gvout_all3,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);
    //    else if(temp ==4) MPI_Reduce(gvout4,gvout_all4,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);
    //    else if(temp ==5) MPI_Reduce(gvout5,gvout_all5,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);
    //    else if(temp ==6) MPI_Reduce(gvout6,gvout_all6,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);

    if (me == 0) {
      fptr = fopen(filename[temp], "w");
      fprintf(fptr, "MDrho from DFT-CES\n");
      fprintf(fptr, "ZYX, Bohr unit\n");
      fprintf(fptr, "% 5d    0.000000    0.000000    0.000000\n", natoms);
      fprintf(fptr, "% 5d% 12.6lf    0.000000    0.000000\n", gnx, gx[0] / 0.52917721);
      fprintf(fptr, "% 5d    0.000000% 12.6lf    0.000000\n", gny, gy[1] / 0.52917721);
      fprintf(fptr, "% 5d    0.000000    0.000000% 12.6lf\n", gnz, gz[2] / 0.52917721);
      for (i = 0; i < natoms; i++) {
        fprintf(fptr, "% 5d% 12.6lf% 12.6lf% 12.6lf% 12.6lf\n", atomns[i], (double)atomns[i], basis[i][0] / 0.52917721, basis[i][1] / 0.52917721, basis[i][2] / 0.52917721);
      }
      for (i = 0; i < gnx; i++) {
        for (j = 0; j < gny; j++) {
          cnt = 0;
          for (k = 0; k < gnz; k++) {
            fprintf(fptr, "% 13.5lE", gvout_all[block * temp + k + j * gnz + i * gnz * gny] * pow(0.52917721, 3) / (nsteps + 1));
            //            if(temp ==0) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            //            else if(temp ==1) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            //            else if(temp ==2) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            //            else if(temp ==3) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            //            else if(temp ==4) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            //            else if(temp ==5) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            //            else if(temp ==6) fprintf(fptr,"% 13.5lE",gvout_all[block*temp+k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
            if (cnt % 6 == 5 && k < gnz - 1)
              fprintf(fptr, "\n");
            cnt++;
          }
          fprintf(fptr, "\n");
        }
      }
      fclose(fptr);
      // printf("\n#######\nDFT-CES: MDrho.cube has been successfully saved\n\n");
      printf("\n#######\nDFT-CES: %s has been successfully saved\n\n", filename[temp]);
    }
  }
}

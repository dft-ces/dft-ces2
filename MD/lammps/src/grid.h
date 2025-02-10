/* -*- c++ -*- ----------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#ifndef LMP_GRID_H
#define LMP_GRID_H

#include "pointers.h"

namespace LAMMPS_NS {

class Grid : protected Pointers {
 public:
  int ncubes;    		                                     // number of output cubefiles
  int nrhocubes;    		                                     // number of input QM rho cubefiles after gaussian convolution
  int savetag;			                                     // save CES grid or not	
  char **savedcube;                                                  // cubefile for save
  char **inputrhocube;                                               // cubefile for save
  int gnx, gny, gnz;                                                 // # of grid points in 3 dim
  double gx[3],gy[3],gz[3];                                          // grid spacing vectors of (unit: Ang)
  int natoms;                                                        // # of atoms
	int tip4p_CES;
  int *atomns;                                                       // atomic numbers for each atom
  double **basis;                                                    // cartesian coords of each atom (unit: Ang)
                                                                     // within unit cell (0 <= coord < 1)
  double *gvin;                                                      // grid values from QM pot (unit: Ry)
//  double *gvout0,*gvout1,*gvout2,*gvout3,*gvout4,*gvout5,*gvout6;                    // grid values for MD rho  (unit: Ry)
//  double *gvout_all0,*gvout_all1,*gvout_all2,*gvout_all3,*gvout_all4,*gvout_all5,*gvout_all6;// grid values for MD rho  (unit: Ry)
  double *gvout,*gvout_all;
//  double *cvin0,*cvin1,*cvin2,*cvin3,*cvin4;                                                      // grid values from QM rho (unit: ebohr-3)
  double *cvin;
  Grid(class LAMMPS *, int, char **);
  ~Grid();
  void read_header(char *);
  void read_content(char *);
  void read_chd(char **);
  void save_grid(char **, int);
};

}

#endif


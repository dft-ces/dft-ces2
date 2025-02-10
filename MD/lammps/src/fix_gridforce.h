/* -*- c++ -*- ----------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gridforce,FixGridForce)

#else

#ifndef LMP_FIX_GRIDFORCE_H
#define LMP_FIX_GRIDFORCE_H

#include "fix.h"
#include <map>
#include <string>

namespace LAMMPS_NS {

class FixGridForce : public Fix {
 public:
  FixGridForce(class LAMMPS *, int, char **);
  ~FixGridForce();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);
  void triInterRep(double, double, double, double, double*, int);
//  void triInter(double, double, double, double, double*, int);
  void triInterInduce(double, double,double, double, double, double*, int);
  double triInterValue(double,double,double,double,double,double,double,double,double,double,double);
  double gridValue(double*, int, int, int);
  void gridGrad(double*, int, int, int, double*);
  void sc_gridGrad(double*, int, int, int, double*);
  int tip4p_CES; //DFT-CES with tip4p tag 
	int find(char *);

//  void grid2ndGrad(double*, int, int, int, double*);

 private:
  double weight;      // weight factor for grid energy and force, usually -1 for QE coupling
  int sfactor;        // supercell factor
  int cubeID;        // cubefile id 
  int cubeID_tip4p_O;        // cubefile id 
  double tmp3[3], tmp4[4],tmp5[4],tmp6[8];
  double newsite[3];
  //double tmp3[3], tmp4[4],tmp5[8];
 
  int force_flag;
//  int force_flag2;

  double Egrid[3];
  double Egrid_all[3];
//  double EgridRep;
//  double Egrid_all, EgridRep_all;
  

 protected:
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  double theta,blen;             // angle and bond types of TIP4P water
  double alpha;                // geometric constraint parameter for TIP4P
  double qdist;
  double RepX;					// Pauli Repulsion parameter of atomic species X in kcal/mol bohr3
  double Polar;					// isotropic polarizability
	char *RepXmeta;

  void compute_newsite(double *, double *, double *, double *);

};

}

#endif
#endif

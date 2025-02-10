#/* ----------------------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
	 Modified by T. Jang
   Copyright (C) 2024 M-design group @ KAIST
------------------------------------------------------------------------- */

#include "fix_gridforce.h"
#include "angle.h"
#include "atom.h"
#include "atom_masks.h"
#include "bond.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "grid.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

/* ---------------------------------------------------------------------- */

FixGridForce::FixGridForce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg != 14 && narg != 7)
    error->all(FLERR, "Illegal fix gridforce command: weight supercellfactor cube-ID element keywords");
  else if (narg == 14 &&  (strcmp(arg[7], "tip4p") != 0) && (strcmp(arg[7], "TIP4P") != 0) )
    error->all(FLERR, "Illegal fix gridforce command: TIP4P arguments - otype htype bdistance angle qdist indicator");
  scalar_flag = 1;
  vector_flag = 1;  
  extvector = 1;    
  size_vector = 3; 
  extscalar = 1;

  weight = atof(arg[3]);
  sfactor = atoi(arg[4]);
  cubeID = atoi(arg[5]);
  RepXmeta = arg[6];

  if (comm->me == 0)
    printf("DFT-CES: weight factor    = %f\n", weight);
  if (comm->me == 0)
    printf("DFT-CES: supercell factor = %d\n", sfactor);
  if (comm->me == 0)
    printf("DFT-CES: cubefile ID      = %d\n", cubeID);

  force_flag = 0;
  Egrid[0] = Egrid[1] = Egrid[2] = 0.0;
  tmp5[0] = tmp5[1] = tmp5[2] = tmp5[3] = 0;
  tmp6[0] = tmp6[1] = tmp6[2] = tmp6[3] = tmp6[4] = tmp6[5] = tmp6[6] = tmp6[7] = 0;
  tip4p_CES = 0;

	if ( narg > 7){
		if ( (strcmp(arg[7], "tip4p") == 0) || (strcmp(arg[7], "TIP4P") == 0) ) {
    tip4p_CES = 1; // DFT-CES with tip4p turn on
    typeO = atoi(arg[8]);
    typeH = atoi(arg[9]);
    blen = atof(arg[10]);
    theta = atof(arg[11]);
    qdist = atof(arg[12]);

    // set alpha parameter
    alpha = qdist / (cos(0.5 * theta * M_PI / 180) * blen);
    cubeID_tip4p_O = atoi(arg[13]);
  } 
  }

	find(RepXmeta);
	
}

/* ---------------------------------------------------------------------- */

FixGridForce::~FixGridForce() {
  return;
}

/* ---------------------------------------------------------------------- */

int FixGridForce::setmask() {
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGridForce::init() {
  // check variables
  if (domain->grid->natoms == 0)
    error->all(FLERR, "Grid data is unavilable");
	if (domain->grid->savetag == 1 && tip4p_CES == 1){
		if (domain->grid->ncubes <= (domain->grid->nrhocubes -1 +3 +1 )) // exclude potentialcube, type of atoms (nrhocubes), dipx-, y-, z- cubes, and tip4p O) 
    	error->all(FLERR, "Grid data is unavilable, allocate more memory to Grid for oxgyen in TIP4P.");
	}

}

/* ---------------------------------------------------------------------- */

void FixGridForce::setup(int vflag) {
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGridForce::post_force(int vflag) {
  int iH1, iH2;
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag; 
  double *q = atom->q;
  int *type = atom->type; // TIP4P
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double fx, fy, fz;
  double repfx, repfy, repfz;

  newsite[0] = newsite[1] = newsite[2] = 0;
  Egrid[0] = Egrid[1] = Egrid[2] = 0.0;
  force_flag = 0;
  tmp5[0] = tmp5[1] = tmp5[2] = tmp5[3] = 0;
  tmp6[0] = tmp6[1] = tmp6[2] = tmp6[3] = tmp6[4] = tmp6[5] = tmp6[6] = tmp6[7] = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (tip4p_CES == 1) {
        if (type[i] == typeO) {
          iH1 = atom->map(tag[i] + 1);
          iH2 = atom->map(tag[i] + 2);
          iH1 = domain->closest_image(i, iH1);
          iH2 = domain->closest_image(i, iH2);
          compute_newsite(x[i], x[iH1], x[iH2], newsite);
          triInterRep(RepX, x[i][0], x[i][1], x[i][2], tmp5, cubeID); // newsite is valid only for the electrostatic force, not Pauli rep and disp
          triInterInduce(q[i], Polar, newsite[0], newsite[1], newsite[2], tmp6, cubeID);
          // save at updated position and force calculate
          Egrid[0] += weight * 23.06092 * q[i] * tmp6[0]; // weight = -1, 1 eV = 23.06092 kcal/mol
          fx = weight * 23.06092 * q[i] * tmp6[1];
          fy = weight * 23.06092 * q[i] * tmp6[2];
          fz = weight * 23.06092 * q[i] * tmp6[3];
          Egrid[1] += 0.237154 * tmp6[4]; // , (0.52917721/13.60569/2)^2 *627.509 = 0.237154
          fx += weight * 0.237154 * tmp6[5];
          fy += weight * 0.237154 * tmp6[6];
          fz += weight * 0.237154 * tmp6[7];
          Egrid[2] += RepX * tmp5[0]; 
          repfx = RepX * tmp5[1];
          repfy = RepX * tmp5[2];
          repfz = RepX * tmp5[3];
          f[i][0] += (repfx + fx) * (1.0 - alpha);
          f[i][1] += (repfy + fy) * (1.0 - alpha);
          f[i][2] += (repfz + fz) * (1.0 - alpha);
          f[iH1][0] += 0.5 * alpha * (fx + repfx);
          f[iH1][1] += 0.5 * alpha * (fy + repfy);
          f[iH1][2] += 0.5 * alpha * (fz + repfz);
          f[iH2][0] += 0.5 * alpha * (fx + repfx);
          f[iH2][1] += 0.5 * alpha * (fy + repfy);
          f[iH2][2] += 0.5 * alpha * (fz + repfz);
          // distribute force at massless point at O,H1,H2
        } else { // type[i]==typeH
          triInterRep(RepX, x[i][0], x[i][1], x[i][2], tmp5, cubeID);
          triInterInduce(q[i], Polar, x[i][0], x[i][1], x[i][2], tmp6, cubeID);
          Egrid[0] += weight * 23.06092 * q[i] * tmp6[0];
          fx = weight * 23.06092 * q[i] * tmp6[1];
          fy = weight * 23.06092 * q[i] * tmp6[2];
          fz = weight * 23.06092 * q[i] * tmp6[3];
          Egrid[1] += 0.237154 * tmp6[4]; 
          fx += weight * 0.237154 * tmp6[5];
          fy += weight * 0.237154 * tmp6[6];
          fz += weight * 0.237154 * tmp6[7];
          Egrid[2] += RepX * tmp5[0];
          repfx = RepX * tmp5[1];
          repfy = RepX * tmp5[2];
          repfz = RepX * tmp5[3];
          f[i][0] += fx + repfx;
          f[i][1] += fy + repfy;
          f[i][2] += fz + repfz;
        }
      } else {
        triInterRep(RepX, x[i][0], x[i][1], x[i][2], tmp5, cubeID);
        triInterInduce(q[i], Polar, x[i][0], x[i][1], x[i][2], tmp6, cubeID);
        Egrid[0] += weight * 23.06092 * q[i] * tmp6[0];
        fx = weight * 23.06092 * q[i] * tmp6[1];
        fy = weight * 23.06092 * q[i] * tmp6[2];
        fz = weight * 23.06092 * q[i] * tmp6[3];
        Egrid[1] += 0.237154 * tmp6[4]; 
        fx += weight * 0.237154 * tmp6[5];
        fy += weight * 0.237154 * tmp6[6];
        fz += weight * 0.237154 * tmp6[7];
        Egrid[2] += RepX * tmp5[0];
        repfx = RepX * tmp5[1];
        repfy = RepX * tmp5[2];
        repfz = RepX * tmp5[3];
        f[i][0] += fx + repfx;
        f[i][1] += fy + repfy;
        f[i][2] += fz + repfz;
      }
    }
  }
}

double FixGridForce::compute_scalar() {
  // only sum across procs one time
  if (force_flag == 0) {
    MPI_Allreduce(&Egrid[0], &Egrid_all[0], 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&Egrid[1], &Egrid_all[1], 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&Egrid[2], &Egrid_all[2], 1, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return Egrid_all[0] + Egrid_all[1] + Egrid_all[2];
}

double FixGridForce::compute_vector(int n) {
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(&Egrid[0], &Egrid_all[0], 3, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return Egrid_all[n]; // an element
}

void FixGridForce::triInterRep(double RepX, double x, double y, double z, double *tmp5, int CUBEID) {
  int i;
  double check;
  double *gvin = domain->grid->gvin;
  double *cvin = NULL;
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  int cube = domain->grid->nrhocubes;
	int savetag = domain->grid->savetag;
  cvin = &(domain->grid->cvin[gnx * gny * gnz * CUBEID]);
  double px, py, pz, xd, yd, zd;
  double pPR000, pPR100, pPR010, pPR001, pPR110, pPR101, pPR011, pPR111;
  double gG[3];
  double ePR000[3], ePR100[3], ePR010[3], ePR001[3], ePR110[3], ePR101[3], ePR011[3], ePR111[3];
  double gsx = domain->grid->gx[0]; // unit in Angstrom;
  double gsy = domain->grid->gy[1];
  double gsz = domain->grid->gz[2];
  double gvol;
  if (x < 0) {
    px = (fmod(x, gnx * gsx) + gnx * gsx) / gsx;
  } else {
    px = fmod(x, gnx * gsx) / gsx;
  }
  if (y < 0) {
    py = (fmod(y, gny * gsy) + gny * gsy) / gsy;
  } else {
    py = fmod(y, gny * gsy) / gsy;
  }
  if (z < 0) {
    pz = (fmod(z, gnz * gsz) + gnz * gsz) / gsz;
  } else {
    pz = fmod(z, gnz * gsz) / gsz;
  }
  xd = (double)(px - (int)px);
  yd = (double)(py - (int)py);
  zd = (double)(pz - (int)pz);

  tmp5[0] = tmp5[1] = tmp5[2] = tmp5[3] = 0;
  pPR000 = gridValue(cvin, (int)px, (int)py, (int)pz); 
  pPR100 = gridValue(cvin, 1 + (int)px, (int)py, (int)pz);
  pPR010 = gridValue(cvin, (int)px, 1 + (int)py, (int)pz);
  pPR001 = gridValue(cvin, (int)px, (int)py, 1 + (int)pz);
  pPR110 = gridValue(cvin, 1 + (int)px, 1 + (int)py, (int)pz);
  pPR101 = gridValue(cvin, 1 + (int)px, (int)py, 1 + (int)pz);
  pPR011 = gridValue(cvin, (int)px, 1 + (int)py, 1 + (int)pz);
  pPR111 = gridValue(cvin, 1 + (int)px, 1 + (int)py, 1 + (int)pz);
  tmp5[0] = triInterValue(pPR000, pPR100, pPR010, pPR001, pPR110, pPR101, pPR011, pPR111, xd, yd, zd);
  check = sqrt(tmp5[0] * tmp5[0]);
  if (check > tmp5[0] || tmp5[0] < 0.000000001) {
    tmp5[0] = 0.0;
  }

  gridGrad(cvin, ((int)px) % gnx, ((int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR000[i] = gG[i];
  gridGrad(cvin, (1 + (int)px) % gnx, ((int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR100[i] = gG[i];
  gridGrad(cvin, ((int)px) % gnx, (1 + (int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR010[i] = gG[i];
  gridGrad(cvin, ((int)px) % gnx, ((int)py) % gny, (1 + (int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR001[i] = gG[i];
  gridGrad(cvin, (1 + (int)px) % gnx, (1 + (int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR110[i] = gG[i];
  gridGrad(cvin, (1 + (int)px) % gnx, ((int)py) % gny, (1 + (int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR101[i] = gG[i];
  gridGrad(cvin, ((int)px) % gnx, (1 + (int)py) % gny, (1 + (int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR011[i] = gG[i];
  gridGrad(cvin, (1 + (int)px) % gnx, (1 + (int)py) % gny, (1 + (int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    ePR111[i] = gG[i];
  tmp5[1] = triInterValue(ePR000[0], ePR100[0], ePR010[0], ePR001[0], ePR110[0], ePR101[0], ePR011[0], ePR111[0], xd, yd, zd);
  tmp5[2] = triInterValue(ePR000[1], ePR100[1], ePR010[1], ePR001[1], ePR110[1], ePR101[1], ePR011[1], ePR111[1], xd, yd, zd);
  tmp5[3] = triInterValue(ePR000[2], ePR100[2], ePR010[2], ePR001[2], ePR110[2], ePR101[2], ePR011[2], ePR111[2], xd, yd, zd);

    if (tip4p_CES == 1 && savetag ==1) {
    if (cubeID == cubeID_tip4p_O){
			int cuber = domain->grid->nrhocubes;
			double *gvout = NULL;
			gvout = &(domain->grid->gvout[gnz * gny * gnx * (cuber + 3)]);
			RepX *= 0.00318720528; // [=] kcal mol-1*e-1*bohr3 to Ry*e-1*bohr3
			pPR000 = (1 - xd) * (1 - yd) * (1 - zd) * abs(weight) * RepX;
			pPR100 = xd * (1 - yd) * (1 - zd) * abs(weight) * RepX;
			pPR010 = (1 - xd) * yd * (1 - zd) * abs(weight) * RepX;
			pPR110 = xd * yd * (1 - zd) * abs(weight) * RepX;
			pPR001 = (1 - xd) * (1 - yd) * zd * abs(weight) * RepX;
			pPR101 = xd * (1 - yd) * zd * abs(weight) * RepX;
			pPR011 = (1 - xd) * yd * zd * abs(weight) * RepX;
		  pPR111 = xd * yd * zd * abs(weight) * RepX;
  		gvol = gsx * gsy * gsz; // Ang3

    	gvout[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pPR000 / gvol / sfactor;
    	gvout[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pPR100 / gvol / sfactor;
    	gvout[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pPR010 / gvol / sfactor;
    	gvout[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pPR001 / gvol / sfactor;
    	gvout[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pPR110 / gvol / sfactor;
    	gvout[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pPR101 / gvol / sfactor;
    	gvout[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pPR011 / gvol / sfactor;
    	gvout[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pPR111 / gvol / sfactor;
		}
  }
  return;
}

void FixGridForce::triInterInduce(double q, double Polar, double x, double y, double z, double *tmp6, int CUBEID) {
  int i;
  double check;
  double dipole[3];
  double *gvin = domain->grid->gvin;
  double px, py, pz, xd, yd, zd;
  double p000, p100, p010, p001, p110, p101, p011, p111;
  double pi000[3], pi100[3], pi010[3], pi001[3], pi110[3], pi101[3], pi011[3], pi111[3]; // induced dipole derivatives of squared field in grid points
  double gG[3];
  double e000[3], e100[3], e010[3], e001[3], e110[3], e101[3], e011[3], e111[3];
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
	int savetag = domain->grid->savetag;
  double gsx = domain->grid->gx[0]; // unit in Angstrom
  double gsy = domain->grid->gy[1];
  double gsz = domain->grid->gz[2];
  double gvol;
  if (x < 0) {
    px = (fmod(x, gnx * gsx) + gnx * gsx) / gsx;
  } else {
    px = fmod(x, gnx * gsx) / gsx;
  }
  if (y < 0) {
    py = (fmod(y, gny * gsy) + gny * gsy) / gsy;
  } else {
    py = fmod(y, gny * gsy) / gsy;
  }
  if (z < 0) {
    pz = (fmod(z, gnz * gsz) + gnz * gsz) / gsz;
  } else {
    pz = fmod(z, gnz * gsz) / gsz;
  }
  xd = (double)(px - (int)px);
  yd = (double)(py - (int)py);
  zd = (double)(pz - (int)pz);
  double Polar2 = Polar / 2; // scaling of alpha to include self-polarization energy

  tmp6[0] = tmp6[1] = tmp6[2] = tmp6[3] = tmp6[4] = tmp6[5] = tmp6[6] = tmp6[7] = 0;

  p000 = gridValue(gvin, (int)px, (int)py, (int)pz);
  p100 = gridValue(gvin, 1 + (int)px, (int)py, (int)pz);
  p010 = gridValue(gvin, (int)px, 1 + (int)py, (int)pz);
  p001 = gridValue(gvin, (int)px, (int)py, 1 + (int)pz);
  p110 = gridValue(gvin, 1 + (int)px, 1 + (int)py, (int)pz);
  p101 = gridValue(gvin, 1 + (int)px, (int)py, 1 + (int)pz);
  p011 = gridValue(gvin, (int)px, 1 + (int)py, 1 + (int)pz);
  p111 = gridValue(gvin, 1 + (int)px, 1 + (int)py, 1 + (int)pz);
  tmp6[0] = triInterValue(p000, p100, p010, p001, p110, p101, p011, p111, xd, yd, zd);
  gridGrad(gvin, ((int)px) % gnx, ((int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    e000[i] = gG[i];
  gridGrad(gvin, (int(1 + px)) % gnx, ((int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    e100[i] = gG[i];
  gridGrad(gvin, ((int)px) % gnx, (int(1 + py)) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    e010[i] = gG[i];
  gridGrad(gvin, ((int)px) % gnx, ((int)py) % gny, (int(1 + pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    e001[i] = gG[i];
  gridGrad(gvin, (int(1 + px)) % gnx, (int(1 + py)) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    e110[i] = gG[i];
  gridGrad(gvin, (int(1 + px)) % gnx, ((int)py) % gny, (int(1 + pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    e101[i] = gG[i];
  gridGrad(gvin, ((int)px) % gnx, (int(1 + py)) % gny, (int(1 + pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    e011[i] = gG[i];
  gridGrad(gvin, (int(1 + px)) % gnx, (int(1 + py)) % gny, (int(1 + pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    e111[i] = gG[i];
  tmp6[1] = triInterValue(e000[0], e100[0], e010[0], e001[0], e110[0], e101[0], e011[0], e111[0], xd, yd, zd);
  tmp6[2] = triInterValue(e000[1], e100[1], e010[1], e001[1], e110[1], e101[1], e011[1], e111[1], xd, yd, zd);
  tmp6[3] = triInterValue(e000[2], e100[2], e010[2], e001[2], e110[2], e101[2], e011[2], e111[2], xd, yd, zd);
  tmp6[4] = -1 * (tmp6[1] * tmp6[1] + tmp6[2] * tmp6[2] + tmp6[3] * tmp6[3]) * Polar2;

  sc_gridGrad(gvin, ((int)px) % gnx, ((int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi000[i] = gG[i];
  sc_gridGrad(gvin, (int(1 + px)) % gnx, ((int)py) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi100[i] = gG[i];
  sc_gridGrad(gvin, ((int)px) % gnx, (int(1+py)) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi010[i] = gG[i];
  sc_gridGrad(gvin, ((int)px) % gnx, ((int)py) % gny, (int(1+pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi001[i] = gG[i];
  sc_gridGrad(gvin, (int(1+px)) % gnx, (int(1+py)) % gny, ((int)pz) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi110[i] = gG[i];
  sc_gridGrad(gvin, (int(1+px)) % gnx, ((int)py) % gny, (int(1+pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi101[i] = gG[i];
  sc_gridGrad(gvin, ((int)px) % gnx, (int(1+py)) % gny, (int(1+pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi011[i] = gG[i];
  sc_gridGrad(gvin, (int(1+px)) % gnx, (int(1+py)) % gny, (int(1+pz)) % gnz, gG);
  for (i = 0; i < 3; i++)
    pi111[i] = gG[i];

	tmp6[5] = Polar2 * (-1) * tmp6[1] * triInterValue(pi000[0], pi100[0], pi010[0], pi001[0], pi110[0], pi101[0], pi011[0], pi111[0], xd, yd, zd); 
  tmp6[6] = Polar2 * (-1) * tmp6[2] * triInterValue(pi000[1], pi100[1], pi010[1], pi001[1], pi110[1], pi101[1], pi011[1], pi111[1], xd, yd, zd);
  tmp6[7] = Polar2 * (-1) * tmp6[3] * triInterValue(pi000[2], pi100[2], pi010[2], pi001[2], pi110[2], pi101[2], pi011[2], pi111[2], xd, yd, zd);

	if (savetag==1){
  // inverse trilinear interpolation for saving MD rho
  p000 = (1 - xd) * (1 - yd) * (1 - zd) * abs(weight) * q;
  p100 = xd * (1 - yd) * (1 - zd) * abs(weight) * q;
  p010 = (1 - xd) * yd * (1 - zd) * abs(weight) * q;
  p110 = xd * yd * (1 - zd) * abs(weight) * q;
  p001 = (1 - xd) * (1 - yd) * zd * abs(weight) * q;
  p101 = xd * (1 - yd) * zd * abs(weight) * q;
  p011 = (1 - xd) * yd * zd * abs(weight) * q;
  p111 = xd * yd * zd * abs(weight) * q;

  dipole[0] = tmp6[1] * Polar2 * abs(weight) * 0.0194469; // 0.52917721/2/13.60569253 = 0.01944690 ebohr
  pi000[0] = (1 - xd) * (1 - yd) * (1 - zd) * dipole[0];
  pi100[0] = xd * (1 - yd) * (1 - zd) * dipole[0];
  pi010[0] = (1 - xd) * yd * (1 - zd) * dipole[0];
  pi110[0] = xd * yd * (1 - zd) * dipole[0];
  pi001[0] = (1 - xd) * (1 - yd) * dipole[0];
  pi101[0] = xd * (1 - yd) * zd * dipole[0];
  pi011[0] = (1 - xd) * yd * zd * dipole[0];
  pi111[0] = xd * yd * zd * dipole[0];

  dipole[1] = tmp6[2] * Polar2 * abs(weight) * 0.01944690;
  pi000[1] = (1 - xd) * (1 - yd) * (1 - zd) * dipole[1];
  pi100[1] = xd * (1 - yd) * (1 - zd) * dipole[1];
  pi010[1] = (1 - xd) * yd * (1 - zd) * dipole[1];
  pi110[1] = xd * yd * (1 - zd) * dipole[1];
  pi001[1] = (1 - xd) * (1 - yd) * dipole[1];
  pi101[1] = xd * (1 - yd) * zd * dipole[1];
  pi011[1] = (1 - xd) * yd * zd * dipole[1];
  pi111[1] = xd * yd * zd * dipole[1];

  dipole[2] = tmp6[3] * Polar2 * abs(weight) * 0.01944690;
  pi000[2] = (1 - xd) * (1 - yd) * (1 - zd) * dipole[2];
  pi100[2] = xd * (1 - yd) * (1 - zd) * dipole[2];
  pi010[2] = (1 - xd) * yd * (1 - zd) * dipole[2];
  pi110[2] = xd * yd * (1 - zd) * dipole[2];
  pi001[2] = (1 - xd) * (1 - yd) * dipole[2];
  pi101[2] = xd * (1 - yd) * zd * dipole[2];
  pi011[2] = (1 - xd) * yd * zd * dipole[2];
  pi111[2] = xd * yd * zd * dipole[2];

  gvol = gsx * gsy * gsz; // Ang3

  double *gvout = NULL;
  gvout = &(domain->grid->gvout[gnz * gny * gnx * CUBEID]);
  gvout[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += p000 / gvol / sfactor;
  gvout[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += p100 / gvol / sfactor;
  gvout[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += p010 / gvol / sfactor;
  gvout[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += p001 / gvol / sfactor;
  gvout[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += p110 / gvol / sfactor;
  gvout[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += p101 / gvol / sfactor;
  gvout[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += p011 / gvol / sfactor;
  gvout[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += p111 / gvol / sfactor;

  int cube = domain->grid->nrhocubes;
  double *gvoutx = NULL;
  double *gvouty = NULL;
  double *gvoutz = NULL;
  gvoutx = &(domain->grid->gvout[gnz * gny * gnx * cube]);
  gvouty = &(domain->grid->gvout[gnz * gny * gnx * (cube + 1)]);
  gvoutz = &(domain->grid->gvout[gnz * gny * gnx * (cube + 2)]);
  gvoutx[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi000[0] / gvol / sfactor;
  gvoutx[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi100[0] / gvol / sfactor;
  gvoutx[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi010[0] / gvol / sfactor;
  gvoutx[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi001[0] / gvol / sfactor;
  gvoutx[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi110[0] / gvol / sfactor;
  gvoutx[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi101[0] / gvol / sfactor;
  gvoutx[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi011[0] / gvol / sfactor;
  gvoutx[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi111[0] / gvol / sfactor;

  gvouty[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi000[1] / gvol / sfactor;
  gvouty[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi100[1] / gvol / sfactor;
  gvouty[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi010[1] / gvol / sfactor;
  gvouty[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi001[1] / gvol / sfactor;
  gvouty[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi110[1] / gvol / sfactor;
  gvouty[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi101[1] / gvol / sfactor;
  gvouty[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi011[1] / gvol / sfactor;
  gvouty[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi111[1] / gvol / sfactor;

  gvoutz[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi000[2] / gvol / sfactor;
  gvoutz[((int)pz) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi100[2] / gvol / sfactor;
  gvoutz[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi010[2] / gvol / sfactor;
  gvoutz[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi001[2] / gvol / sfactor;
  gvoutz[((int)pz) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi110[2] / gvol / sfactor;
  gvoutz[((int)pz + 1) % gnz + (((int)py) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi101[2] / gvol / sfactor;
  gvoutz[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px) % gnx) * gnz * gny] += pi011[2] / gvol / sfactor;
  gvoutz[((int)pz + 1) % gnz + (((int)py + 1) % gny) * gnz + (((int)px + 1) % gnx) * gnz * gny] += pi111[2] / gvol / sfactor;

	}
  return;
}

double FixGridForce::triInterValue(double c000, double c100, double c010, double c001, double c110, double c101, double c011, double c111, double xd, double yd, double zd) {
  double c00, c10, c01, c11, c0, c1, c;
  c00 = c000 * (1 - xd) + c100 * xd;
  c01 = c001 * (1 - xd) + c101 * xd;
  c10 = c010 * (1 - xd) + c110 * xd;
  c11 = c011 * (1 - xd) + c111 * xd;
  c0 = c00 * (1 - yd) + c10 * yd;
  c1 = c01 * (1 - yd) + c11 * yd;
  c = c0 * (1 - zd) + c1 * zd;
  return c;
}

double FixGridForce::gridValue(double *gvin, int x, int y, int z) {
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  double gValue = 0;

  gValue = gvin[(z % gnz) + (y % gny) * gnz + (x % gnx) * gnz * gny];

  return gValue;
}

void FixGridForce::sc_gridGrad(double *gvin, int x, int y, int z, double *gG) { // gradient of squared field output
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  double gsx = domain->grid->gx[0];
  double gsy = domain->grid->gy[1];
  double gsz = domain->grid->gz[2];
  double xbef, xaft, ybef, yaft, zbef, zaft;
  int xplus, xminus, yplus, yminus, zplus, zminus;
  double nzz[3], pzz[3], znz[3], zpz[3], zzn[3], zzp[3];
  double ingG[3];
  ingG[0] = ingG[1] = ingG[2] = 0;
  gG[0] = gG[1] = gG[2] = 0;
	xplus = (x+1)%gnx;
	xminus = (gnx + x -1)%gnx;
	yplus = (y+1)%gny;
	yminus = (gny + y -1)%gny;
	zplus = (z+1)%gnz;
	zminus = (gnz + z -1)%gnz;
  // c-100=nzz
  xbef = gvin[z + y * gnz + ((gnx + x - 2) % gnx) * gnz * gny];
  xaft = gvin[z + y * gnz + x * gnz * gny];
  ybef = gvin[z + ((gny + y - 1) % gny) * gnz + xminus * gnz * gny];
  yaft = gvin[z + ((y + 1) % gny) * gnz + xminus * gnz * gny];
  zbef = gvin[((gnz + z - 1) % gnz) + y * gnz + xminus * gnz * gny];
  zaft = gvin[(z + 1) % gnz + y * gnz + xminus * gnz * gny];
  nzz[0] = (xaft - xbef) / (2 * gsx);
  nzz[1] = (yaft - ybef) / (2 * gsy);
  nzz[2] = (zaft - zbef) / (2 * gsz);
  // c100=pzz
  xbef = gvin[z + y * gnz + x * gnz * gny];
  xaft = gvin[z + y * gnz + ((x + 2) % gnx) * gnz * gny];
  ybef = gvin[z + ((gny + y - 1) % gny) * gnz + xplus * gnz * gny];
  yaft = gvin[z + ((y + 1) % gny) * gnz + xplus * gnz * gny];
  zbef = gvin[((gnz + z - 1) % gnz) + y * gnz + xplus * gnz * gny];
  zaft = gvin[(z + 1) % gnz + y * gnz + xplus * gnz * gny];
  pzz[0] = (xaft - xbef) / (2 * gsx);
  pzz[1] = (yaft - ybef) / (2 * gsy);
  pzz[2] = (zaft - zbef) / (2 * gsz);
  // c0-10=znz
  xbef = gvin[z + yminus * gnz + ((gnx + x - 1) % gnx) * gnz * gny];
  xaft = gvin[z + yminus * gnz + ((x + 1) % gnx) * gnz * gny];
  ybef = gvin[z + ((gny + y - 2) % gny) * gnz + x * gnz * gny];
  yaft = gvin[z + y * gnz + x * gnz * gny];
  zbef = gvin[((gnz + z - 1) % gnz) + yminus * gnz + x * gnz * gny];
  zaft = gvin[(z + 1) % gnz + yminus * gnz + x * gnz * gny];
  znz[0] = (xaft - xbef) / (2 * gsx);
  znz[1] = (yaft - ybef) / (2 * gsy);
  znz[2] = (zaft - zbef) / (2 * gsz);
  // c010=zpz
  xbef = gvin[z + yplus * gnz + ((gnx + x - 1) % gnx) * gnz * gny];
  xaft = gvin[z + yplus * gnz + ((x + 1) % gnx) * gnz * gny];
  ybef = gvin[z + y * gnz + x * gnz * gny];
  yaft = gvin[z + ((y + 2) % gny) * gnz + x * gnz * gny];
  zbef = gvin[((gnz + z - 1) % gnz) + yplus * gnz + x * gnz * gny];
  zaft = gvin[(z + 1) % gnz + yplus * gnz + x * gnz * gny];
  zpz[0] = (xaft - xbef) / (2 * gsx);
  zpz[1] = (yaft - ybef) / (2 * gsy);
  zpz[2] = (zaft - zbef) / (2 * gsz);
  // c00-1=zzn
  xbef = gvin[zminus + y * gnz + ((gnx + x - 1) % gnx) * gnz * gny];
  xaft = gvin[zminus + y * gnz + ((x + 1) % gnx) * gnz * gny];
  ybef = gvin[zminus + ((gny + y - 1) % gny) * gnz + x * gnz * gny];
  yaft = gvin[zminus + ((y + 1) % gny) * gnz + x * gnz * gny];
  zbef = gvin[((gnz + z - 2) % gnz) + y * gnz + x * gnz * gny];
  zaft = gvin[z + y * gnz + x * gnz * gny];
  zzn[0] = (xaft - xbef) / (2 * gsx);
  zzn[1] = (yaft - ybef) / (2 * gsy);
  zzn[2] = (zaft - zbef) / (2 * gsz);
  // c001=zzp
  xbef = gvin[zplus + y * gnz + ((gnx + x - 1) % gnx) * gnz * gny];
  xaft = gvin[zplus + y * gnz + ((x + 1) % gnx) * gnz * gny];
  ybef = gvin[zplus + ((gny + y - 1) % gny) * gnz + x * gnz * gny];
  yaft = gvin[zplus + ((y + 1) % gny) * gnz + x * gnz * gny];
  zbef = gvin[z + y * gnz + x * gnz * gny];
  zaft = gvin[(z + 2) % gnz + y * gnz + x * gnz * gny];
  zzp[0] = (xaft - xbef) / (2 * gsx);
  zzp[1] = (yaft - ybef) / (2 * gsy);
  zzp[2] = (zaft - zbef) / (2 * gsz);
  gG[0] = -1 * (pzz[0] - nzz[0]) / (2 * gsx);
  gG[1] = -1 * (zpz[1] - znz[1]) / (2 * gsy);
  gG[2] = -1 * (zzp[2] - zzn[2]) / (2 * gsz);
  return;
}

void FixGridForce::gridGrad(double *gvin, int x, int y, int z, double *gG) {
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  double gsx = domain->grid->gx[0];
  double gsy = domain->grid->gy[1];
  double gsz = domain->grid->gz[2];

  gG[0] = gG[1] = gG[2] = 0;

  gG[0] = -1 * (-1 * gvin[z + y * gnz + ((x + 2) % gnx) * gnz * gny] + 8 * gvin[z + y * gnz + ((x + 1) % gnx) * gnz * gny] - 8 * gvin[z + y * gnz + (x ? (x - 1) : (gnx - 1)) * gnz * gny] + gvin[z + y * gnz + ((gnx + x - 2) % gnx) * gnz * gny]) / (12 * gsx);
  gG[1] = -1 * (-1 * gvin[z + ((y + 2) % gny) * gnz + x * gnz * gny] + 8 * gvin[z + ((y + 1) % gny) * gnz + x * gnz * gny] - 8 * gvin[z + (y ? (y - 1) : (gny - 1)) * gnz + x * gnz * gny] + gvin[z + ((gny + y - 2) % gny) * gnz + x * gnz * gny]) / (12 * gsy);
  gG[2] = -1 * (-1 * gvin[((z + 2) % gnz) + y * gnz + x * gnz * gny] + 8 * gvin[((z + 1) % gnz) + y * gnz + x * gnz * gny] - 8 * gvin[(z ? (z - 1) : (gnz - 1)) + y * gnz + x * gnz * gny] + gvin[((gnz + z - 2) % gnz) + y * gnz + x * gnz * gny]) / (12 * gsz);
  // gG[0] = -1*(gvin[z+y*gnz+((x+1)%gnx)*gnz*gny]-gvin[z+y*gnz+(x?(x-1):(gnx-1))*gnz*gny])/(2*gsx);
  // gG[1] = -1*(gvin[z+((y+1)%gny)*gnz+x*gnz*gny]-gvin[z+(y?(y-1):(gny-1))*gnz+x*gnz*gny])/(2*gsy);
  // gG[2] = -1*(gvin[((z+1)%gnz)+y*gnz+x*gnz*gny]-gvin[(z?(z-1):(gnz-1))+y*gnz+x*gnz*gny])/(2*gsz);

  return;
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void FixGridForce::compute_newsite(double *xO, double *xH1,
                                   double *xH2, double *xM) {
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];

  xM[0] = xO[0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = xO[1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = xO[2] + alpha * 0.5 * (delz1 + delz2);
}

int FixGridForce::find(char *name) {
  if(name==NULL) return -1;
  if(strcmp(name,"H") == 0) {RepX=866.6; Polar=1.693;} // unit: real for RepX
  if(strcmp(name,"O") == 0) {RepX=9766.9; Polar=5.20;}
  if(strcmp(name,"N") == 0) {RepX=6593.4; Polar=7.25;}
  if(strcmp(name,"S") == 0) {RepX=2127.6; Polar=19.5;}
  if(strcmp(name,"C") == 0) {RepX=3231.7; Polar=11.7;}
  if(strcmp(name,"Li") == 0) {RepX=7091.; Polar=0.193;}
  if(strcmp(name,"Na") == 0) {RepX=15575.; Polar=0.93;}
  if(strcmp(name,"K") == 0) {RepX=14753.; Polar=5.05;}
  if(strcmp(name,"Rb") == 0) {RepX=12086.; Polar=8.32;}
  if(strcmp(name,"Cs") == 0) {RepX=7047.; Polar=15.0;}
  if(strcmp(name,"He") == 0) {RepX=14922.; Polar=1.38;}
  if(strcmp(name,"Ne") == 0) {RepX=13811.; Polar=2.67;}
  if(strcmp(name,"Ar") == 0) {RepX=6953.; Polar=11.1;}
  if(strcmp(name,"Kr") == 0) {RepX=4780.; Polar=16.8;}
  if(strcmp(name,"Xe") == 0) {RepX=2844.; Polar=27.2;}
  if(strcmp(name,"F") == 0) {RepX=1632.; Polar=15.0;}
  if(strcmp(name,"Cl") == 0) {RepX=1691.; Polar=30.3;}
  if(strcmp(name,"Br") == 0) {RepX=1365.; Polar=42.8;}
  if(strcmp(name,"I") == 0) {RepX=1216.; Polar=61.7;}
  if(strcmp(name,"P") == 0) {RepX=2744.5; Polar=24.8;}
  return -1;
}


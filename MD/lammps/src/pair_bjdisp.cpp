/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eduardo Bringa (LLNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pair_bjdisp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;
#define newb		2.1
#define PI		  3.141592

/* ---------------------------------------------------------------------- */

PairBjdisp::PairBjdisp(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairBjdisp::~PairBjdisp()
{
  if (!copymode) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(a);
    memory->destroy(rho);
    memory->destroy(c);
    memory->destroy(rhoinv);
    memory->destroy(buck1);
    memory->destroy(buck2);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairBjdisp::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair,tmpx,tmpy,tmpz;
  double rsq,r2inv,r6,r6inv,forcebuck,factor_lj;
  double r;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nsteps = update->nsteps;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
				r6 = rsq*rsq*rsq;
        r6inv = 1.0/(r6+rhoinv[itype][jtype]);
        forcebuck = -1*buck2[itype][jtype]*r6inv*r6inv*r6;
			} else forcebuck = 0.0;

      fpair = (factor_lj*forcebuck) * r2inv;
			tmpx = delx*fpair; 
			tmpy = dely*fpair; 
			tmpz = delz*fpair;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      if (newton_pair || j < nlocal) {
				tmpx = delx*fpair;
				tmpy = dely*fpair;
				tmpz = delz*fpair;
				f[j][0] -= tmpx;
				f[j][1] -= tmpy;
				f[j][2] -= tmpz;
			}
      if (eflag) {
        if (rsq < cutsq[itype][jtype]) {
          evdwl = -1*c[itype][jtype]*(1.0/(r6+rhoinv[itype][jtype])) - offset[itype][jtype];
          evdwl *= factor_lj;
        } else evdwl = 0.0;
      }
			ecoul = 0.0;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,ecoul,fpair,delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBjdisp::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(a,n+1,n+1,"pair:a");
  memory->create(rho,n+1,n+1,"pair:rho");
  memory->create(c,n+1,n+1,"pair:c");
  memory->create(rhoinv,n+1,n+1,"pair:rhoinv");
  memory->create(buck1,n+1,n+1,"pair:buck1");
  memory->create(buck2,n+1,n+1,"pair:buck2");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBjdisp::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBjdisp::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double a_one = force->numeric(FLERR,arg[2]);
  double rho_one = force->numeric(FLERR,arg[3]);
  if (rho_one < 0) error->all(FLERR,"Incorrect args for pair coefficients");
  double c_one = force->numeric(FLERR,arg[4]);

  double cut_lj_one = cut_lj_global;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      c[i][j] = c_one*rho_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBjdisp::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBjdisp::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = cut_lj[i][j];
  double unit2angs = 0.52917721; 
  double unit2kcalmol = 13.77929 ; // Eh a0^6 to kcal/mol Angs^6
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  double r0ij = pow((sqrt(2/PI)*a[i][j]/3),1./3);
  rhoinv[i][j] = pow((newb*r0ij),6.0);
  buck1[i][j] = 0;
  buck2[i][j] = 6.0*c[i][j];

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    offset[i][j] = -1*c[i][j]/(pow(cut_lj[i][j],6.0)+rhoinv[i][j]);
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  rho[j][i] = rho[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
	}

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBjdisp::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&rho[i][j],sizeof(double),1,fp);
        fwrite(&c[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBjdisp::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&a[i][j],sizeof(double),1,fp);
          fread(&rho[i][j],sizeof(double),1,fp);
          fread(&c[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rho[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBjdisp::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBjdisp::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBjdisp::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,a[i][i],rho[i][i],c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBjdisp::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,
              a[i][j],rho[i][j],c[i][j],cut_lj[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairBjdisp::single(int i, int j, int itype, int jtype,
                               double rsq,
                               double factor_lj,
                               double &fforce)
{
  double r2inv,r6,r6inv,r,forcebuck,phibuck;

  r2inv = 1.0/rsq;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6 = rsq*rsq*rsq;
    r6inv = 1.0/(r6 + rhoinv[itype][jtype]);
    r = sqrt(rsq);
    forcebuck = -1*buck2[itype][jtype]*r6inv*r6inv*r6;
  } else forcebuck = 0.0;
  fforce = (factor_lj*forcebuck) * r2inv;

  double eng = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    phibuck = -1*c[itype][jtype]*(1.0/(r6+rhoinv[itype][jtype])) -
      offset[itype][jtype];
    eng += factor_lj*phibuck;
  }
  return eng;
}

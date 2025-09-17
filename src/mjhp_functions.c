#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"

void find_MJregion(VectorIndices *VECT, UnitCell *UC) 
{
  /*This function finds the shell in reciprocal space that both planewaves
 * of a pair must lie within to contribute to the mjhp*/
  int H, K, L; /*HKL indices*/
  double x0, y0, z0; /*mid points of HKL indices*/
  double kx, ky, kz; /*cartesian coordinates of midpoints*/
  double mag_Ghkl; /*magnitude of midpoint*/
  double rad; /*value to add and subtract from midpoint magnitude to define a shell*/
  double minr, maxr; /*smaller and larger sphere that define the shell in rec space*/
  
  printf( "\nFinding mjhp shell in reciprocal space.\n");

  /*Initialize Variables*/
  H = VECT->H;
  K = VECT->K;
  L = VECT->L;
  rad = 0.075;

  /*Find JZ face-center*/
  x0 = (double) H/2.0;
  y0 = (double) K/2.0;
  z0 = (double) L/2.0;
  /*calculate ang-1 coords of jzfc*/
  kx = x0*UC->ang_ax_star + y0*UC->ang_bx_star + z0*UC->ang_cx_star;
  ky = x0*UC->ang_ay_star + y0*UC->ang_by_star + z0*UC->ang_cy_star;
  kz = x0*UC->ang_az_star + y0*UC->ang_bz_star + z0*UC->ang_cz_star;

  /*calc magnitude/ distance of kx,ky,kx*/
  mag_Ghkl = sqrt(kx*kx+ky*ky+kz*kz);
  printf( "\tCenter of JZ face = %lf %lf %lf\t (%lf %lf %lf ang-1) |G_hkl| = %lf ang-1\n", x0, y0, z0, kx, ky, kz, mag_Ghkl);
  
  minr = mag_Ghkl - rad;
  maxr = mag_Ghkl + rad;
  printf( "\tShell: inner radius = %lf ang-1\t outer radius = %lfang-1\n", minr, maxr);
  VECT->minr = minr;
  VECT->maxr = maxr;

} //END of find_MJregion 

void concatinate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP) 
{
  /*This function combines the local and nonlocal energy to find the total 
 * mjhp potential energy*/
  int nEstep; /*number of energy steps*/
  int dE; /*energy step index*/
  double total_potential; /*total summed potential energy*/

  /*initialize variables*/
  nEstep = ESTP->nEstep;

  /*Allocate Memory*/
  ECON->total_potential = AllocateMemory_oneD_double(ECON->total_potential, nEstep);
  for (dE=0;dE<nEstep;dE++) ECON->total_potential[dE] = 0.0;
  
  printf( "\nConcatinating local and nonlocal potential energy grids.\n");
  /*add local and nonlocal grids to find total potential energy*/
  total_potential = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    ECON->total_potential[dE] = ECON->local[dE] + ECON->nonlocal[dE] - ECON->KE_correction[dE];
    total_potential += ECON->total_potential[dE];
  }
  printf("\tTotal Potential Energy = %lf \n", total_potential);

  /*free memory for local and nonlocal*/
  ECON->local = FreeMemory_oneD_double(ECON->local);
  ECON->nonlocal = FreeMemory_oneD_double(ECON->nonlocal);
  ECON->KE_correction = FreeMemory_oneD_double(ECON->KE_correction);
}

void integrate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP, AtomicVariables * ATM)
{
  /*This function integrates the potential energy up to the Ef*/
  int nEstep; /*number of energy steps*/
  int dE; /*energy step index*/
  int dE_zero; /*zero energy in terms of energy steps*/
  double Eint_Ha; /*integrated energy in Hartrees*/
  double Eint_eV_peratom; /*integrated energy in eV per atom*/
  int natom; /*number of atoms*/
 
  /*initializing variables*/
  nEstep = ESTP->nEstep;
  dE_zero = ESTP->dE_zero;
  natom = ATM->natom;

  printf("\nIntegrating the total potential energy.\n");

  Eint_Ha = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    /*check if energy step is below zero, E_F*/
    if (dE > dE_zero) continue;
    /*if so, store contribution*/
    Eint_Ha += ECON->total_potential[dE];
  }
  /*convert energy to eV and divide by number of atoms*/
  Eint_eV_peratom = (Eint_Ha*HATOEV)/((double)natom);

  printf( "\t iMJHP = %lf eV/atom\n", Eint_eV_peratom);

} //END of integrate_HKL_potential

void integrate_Erange_potential(EnergyContribution * ECON, EnergyStep * ESTP, VectorIndices * VECT)
{
  int nEstep;
  int dE;
  double eval;
  double Eint_Ha;
  double Eint_eV;
  double Elow;
  double scanE_start;  
  double scanE_stop;
  
  nEstep = ESTP->nEstep;
  ESTP->nEstep = nEstep;
  Elow = ESTP->bandE_min*EMESH;
  scanE_start = VECT->scanE_start;
  scanE_stop = VECT->scanE_stop;

  printf("\nIntegrating MJHP from %lf to %lf\n", scanE_start, scanE_stop);
  
  Eint_Ha = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    eval = ((double)dE+Elow)/EMESH;
    if ((eval>scanE_stop)||(eval<scanE_start)) continue;
    Eint_Ha += ECON->local[dE];
  }
  Eint_eV = (Eint_Ha*HATOEV);
  printf("\t iMJHP = %lf eV\n", Eint_eV);

} //END of integrate E range



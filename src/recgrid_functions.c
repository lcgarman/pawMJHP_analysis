#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"

void find_HKLgrid_max(NumberGrid *GRD) 
{  
  /*Calculates the max miller indices based of ngfft grid*/

  int ngfftx, ngffty, ngfftz; /*number grid*/

  /*initialize variables*/
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;

  /*if the nggft is even the max value is half;
 * if it is odd, the max value is the half of nggft-1*/
  if (ngfftx%2==0) {
	GRD->hmax = ngfftx/2;
  } else {
	GRD->hmax = (ngfftx-1)/2;
  }
  if (ngffty%2==0) {
    GRD->kmax = ngffty/2;
  } else {
	GRD->kmax = (ngffty-1)/2;
  }
  if (ngfftz%2==0) {
    GRD->lmax = ngfftz/2;
  } else {
	GRD->lmax = (ngfftz-1)/2;
  }

}  //END of find_HKLgrid_max function

void prepare_HKLgrid(NumberGrid *GRD) 
{
  /*This function unwraps the nggfft grid in reciprocal space 
 * so it is centered about 000*/
  int h0, k0, l0; /*original, positive hkl indices of reciprocal space grid*/
  int h1, k1, l1;  /*true hkl indices of reciprocal space grid*/
  int ngfftx, ngffty, ngfftz; /*number grid*/
  int hmax, kmax, lmax; /*max points on ngfft grids*/

  /*initialize variables*/
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  hmax = GRD->hmax;
  kmax = GRD->kmax;
  lmax = GRD->lmax;

  /*allocate memory for HKL arrays*/
  GRD->h_grid = AllocateMemory_threeD_int(GRD->h_grid, ngfftx, ngffty, ngfftz);
  GRD->k_grid = AllocateMemory_threeD_int(GRD->k_grid, ngfftx, ngffty, ngfftz);
  GRD->l_grid = AllocateMemory_threeD_int(GRD->l_grid, ngfftx, ngffty, ngfftz);

  /*unwrap the ngfft grids to their negative values, so ngfft is centered around 000*/
  for (h0=0; h0<ngfftx; h0++) {
    for (k0=0; k0<ngffty; k0++) {
      for (l0=0; l0<ngfftz; l0++) { 

        if (h0 > hmax) {
          h1 = h0 - ngfftx;
        }  else {
             h1 = h0;
        }
        if (k0 > kmax) {
          k1 = k0 - ngffty;
        }  else {
             k1 = k0;
        }
        if (l0 > lmax) {
          l1 = l0 - ngfftz;
        }  else {
             l1 = l0;
        }
        /*store the grids in grid structure*/
        GRD->h_grid[h0][k0][l0] = h1;
        GRD->k_grid[h0][k0][l0] = k1;
        GRD->l_grid[h0][k0][l0] = l1;
      }
    }
  }

} //END prepare_HKLgrid function


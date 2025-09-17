#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "twotheta_functions.h"
#include "hkl_functions.h"

void calculate_powder_pattern(TwoTheta * TTH, BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC, Symmetry* SYM) 
{
  /*This function calculates a simulated powder diffraction pattern based off of the abinit DEN file*/
  int ngfftx, ngffty, ngfftz; /*number grid*/
  int h0, k0, l0; /*indices for ngfft grid*/
  int h1, k1, l1; /*indices for true reciprocal space grid*/
  double F_hkl_val; /*structure factor*/
  double kx, ky, kz; /*coordinates in inverse ang of reflection*/
  double mag_G; /*magnitude of reflection*/
  double d_hkl_val; /*distance in rec space of reflection*/
  double lambda=0.71073; /*wavelenght of Mo radiation in ang*/
  double rad_theta; /*theta of reflection in radians*/
  double deg_theta; /*theta of reflection in degrees*/
  double twotheta_val; /*two theta of reflection in degrees*/
  int nsym; /*number of symmetry elements*/
  int sym; /*symmetry index*/
  int match; /*integer tracking if a reflection match has been found*/
  int j; /*search index*/
  int symm_nhkl; /*number of symmetry-independent reflections*/
  double F_diff; /*difference between two structure factors of symmetry-related reflections*/
  int ngfft_size; /*size of ngfft grid*/
  int * rflc_mult; /*multiplicity of reflections*/
  int * H_arr; /*array of HKL reflection indices*/
  int * K_arr;
  int * L_arr;
  int H_symm, K_symm, L_symm; /*symmetry-related HKL indices*/
  double tol; /*tolerance on structure factors for related reflections*/
  double * two_theta; /*array of two theta of reflections*/
  double * d_hkl; /*array of distance of reflections in rec space*/
  double * F_hkl; /*array of structure factors*/
  int* hkl_sum; /*sum of HKL indices to organize reflection storing*/
  int hkl_sum_new; 

  /*initialize variables*/
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  nsym = SYM->nsym;
  symm_nhkl = 0;
  tol = 1e-10;
  
  /*nullify local arrays*/
  rflc_mult = NULL; 
  H_arr = NULL; 
  K_arr = NULL;
  L_arr = NULL;
  hkl_sum = NULL;
  two_theta = NULL;
  d_hkl = NULL;
  F_hkl = NULL;
  /*allocate memory for local arrays*/
  ngfft_size = (ngfftx*ngffty*ngfftz);
  rflc_mult = AllocateMemory_oneD_int(rflc_mult, ngfft_size);
  H_arr = AllocateMemory_oneD_int(H_arr, ngfft_size);
  K_arr = AllocateMemory_oneD_int(K_arr, ngfft_size);
  L_arr = AllocateMemory_oneD_int(L_arr, ngfft_size);
  hkl_sum = AllocateMemory_oneD_int(hkl_sum, ngfft_size);
  two_theta = AllocateMemory_oneD_double(two_theta, ngfft_size);
  d_hkl = AllocateMemory_oneD_double(d_hkl, ngfft_size);
  F_hkl = AllocateMemory_oneD_double(F_hkl, ngfft_size);
  /*end of allocation*/

  printf("\nCalculating powder pattern.\n");
  /*Begin finding powder pattern*/
  for (h0=0;h0<ngfftx;h0++) {
	for (k0=0;k0<ngffty;k0++) {
	  for (l0=0;l0<ngfftz;l0++) {

		/*skip the 000 point in reciprocal space*/
		if ((h0==0)&&(k0==0)&&(l0==0)) continue;

		/*find the structure factor for this point*/
		F_hkl_val = gsl_complex_abs2(BIN->rec_grid[h0][k0][l0]);

		/*if the structure factor is zero skip*/
		if (F_hkl_val == 0.0) continue;

        /*assign true HKL reflection values*/
		h1 = GRD->h_grid[h0][k0][l0];
		k1 = GRD->k_grid[h0][k0][l0];
		l1 = GRD->l_grid[h0][k0][l0];

        /*find the coordinates of each reflection in inverse ang*/
		kx = h1*UC->ang_ax_star+k1*UC->ang_bx_star+l1*UC->ang_cx_star;
		ky = h1*UC->ang_ay_star+k1*UC->ang_by_star+l1*UC->ang_cy_star;
		kz = h1*UC->ang_az_star+k1*UC->ang_bz_star+l1*UC->ang_cz_star;

        /*find magnitude of reflection*/
		mag_G = sqrt(kx*kx + ky*ky + kz*kz);
        /*find distance in rec space of the reflection*/
		d_hkl_val = (2.0*PI)/mag_G;
        /*if reflection is too high resolution for this wavelength continue*/
		if (d_hkl_val < 0.5*lambda) continue;
        /*find theta value in radians using braggs law*/ 
		rad_theta = asin((lambda)/(2*d_hkl_val));
        /*convert theta value to degrees*/
		deg_theta = rad_theta * 180.0/PI;
		twotheta_val = deg_theta*2.0;
        /*if twotheta is greater than 60 skip*/
		if (twotheta_val > 60) continue;

        /*find symmetry-related hkl reflections*/
        match = 0;
		printf("\tHKL = %d %d %d\n", h1, k1, l1);
        /*apply symmetry to hkl indices, and find matches to other hkl indices already found*/
        for(sym=0;sym<nsym;sym++) {
          /*if a match was already found discontinue search*/
          if (match == 1) break;
	      H_symm = SYM->symrel[0][0][sym]*h1 + SYM->symrel[1][0][sym]*k1 + SYM->symrel[2][0][sym]*l1;
	      K_symm = SYM->symrel[0][1][sym]*h1 + SYM->symrel[1][1][sym]*k1 + SYM->symrel[2][1][sym]*l1;
    	  L_symm = SYM->symrel[0][2][sym]*h1 + SYM->symrel[1][2][sym]*k1 + SYM->symrel[2][2][sym]*l1;
	      for (j=0;j<symm_nhkl;j++) {
		    F_diff = F_hkl_val - F_hkl[j];
            /*check if the symmetry-related HKL indices match that of a previously found reflection with matching Fhkl values*/ 
            if ((H_symm==H_arr[j])&&(K_symm==K_arr[j])&&(L_symm==L_arr[j])&&(fabs(F_diff)<tol)) {
              /*if match found increase multiplicity of reflection*/
              rflc_mult[j]++;
              match = 1;
              printf("\t\tMATCH(%e): %d %d %d F_hkl=%e\t == %d %d %d F_hkl=%e\n", F_diff, h1, k1, l1, F_hkl_val, H_arr[j], K_arr[j], L_arr[j], F_hkl[j]);
              /*find highest combination of hkl indices*/
              hkl_sum_new = h1+k1+l1;
              if (hkl_sum_new>hkl_sum[j]) {
                H_arr[j] = h1;
                K_arr[j] = k1;
                L_arr[j] = l1;
                hkl_sum[j] = hkl_sum_new;
			  }
              /*if a match was found, break the search for a match*/
              break;
            }
          }
        }
        /*if no match found, store this HKL as new reflection*/
        if (match == 0) {
	      H_arr[symm_nhkl] = h1;
	      K_arr[symm_nhkl] = k1;
		  L_arr[symm_nhkl] = l1;
		  two_theta[symm_nhkl] = twotheta_val;
		  d_hkl[symm_nhkl] = d_hkl_val;
		  F_hkl[symm_nhkl] = F_hkl_val;
		  rflc_mult[symm_nhkl] = 1;
          hkl_sum[symm_nhkl] = h1+k1+l1;
		  symm_nhkl++;
		}
	  }  //END l loop
	}  //END k loop
  }  //END h loop
  
  /*allocate memory for cpp arrays*/
  TTH->nrflc = symm_nhkl;
  TTH->rflc_mult = AllocateMemory_oneD_int(TTH->rflc_mult, TTH->nrflc);
  TTH->rflc_H = AllocateMemory_oneD_int(TTH->rflc_H, TTH->nrflc);
  TTH->rflc_K = AllocateMemory_oneD_int(TTH->rflc_K, TTH->nrflc);
  TTH->rflc_L = AllocateMemory_oneD_int(TTH->rflc_L, TTH->nrflc);
  TTH->two_theta = AllocateMemory_oneD_double(TTH->two_theta, TTH->nrflc);
  TTH->d_hkl = AllocateMemory_oneD_double(TTH->d_hkl, TTH->nrflc);
  TTH->F_hkl = AllocateMemory_oneD_double(TTH->F_hkl, TTH->nrflc);

  printf("\nStoring Reflections:\n");
  /*store vaiables in TTH struct*/
  for (j=0;j<symm_nhkl;j++) {
    TTH->rflc_H[j] = H_arr[j];
    TTH->rflc_K[j] = K_arr[j];
    TTH->rflc_L[j] = L_arr[j];
    TTH->two_theta[j] = two_theta[j];
    TTH->d_hkl[j] = d_hkl[j];
    TTH->F_hkl[j] = F_hkl[j];
    TTH->rflc_mult[j] = rflc_mult[j];
    printf("\t%d HKL = %d %d %d\t mult = %d\t F_hkl = %lf\n", j, H_arr[j], K_arr[j], L_arr[j], rflc_mult[j], F_hkl[j]); 
  }

  /*free local allocated variables*/
  rflc_mult = FreeMemory_oneD_int(rflc_mult);
  H_arr = FreeMemory_oneD_int(H_arr);
  K_arr = FreeMemory_oneD_int(K_arr);
  L_arr = FreeMemory_oneD_int(L_arr);
  hkl_sum = FreeMemory_oneD_int(hkl_sum);
  two_theta = FreeMemory_oneD_double(two_theta);
  d_hkl = FreeMemory_oneD_double(d_hkl);
  F_hkl = FreeMemory_oneD_double(F_hkl);

  /*Free HKL grids dont need again*/
  GRD->h_grid = FreeMemory_threeD_int(GRD->h_grid, ngfftx, ngffty);
  GRD->k_grid = FreeMemory_threeD_int(GRD->k_grid, ngfftx, ngffty);
  GRD->l_grid = FreeMemory_threeD_int(GRD->l_grid, ngfftx, ngffty);
  
} //END of calculate_powder_pattern

void fold_reflections_toBZ(TwoTheta * TTH) 
{
  /*This function wraps the midpoints of the reflections into the 1st Brillouin zone*/
  int nrflc; /*number of reflections*/
  int n; /*reflection index*/

  /*initialize variables*/
  nrflc = TTH->nrflc;

  /*allocate memory for variables*/
  TTH->hpw = AllocateMemory_oneD_int(TTH->hpw, nrflc);
  TTH->kpw = AllocateMemory_oneD_int(TTH->kpw, nrflc);
  TTH->lpw = AllocateMemory_oneD_int(TTH->lpw, nrflc);
  TTH->BZkpt = AllocateMemory_twoD_double(TTH->BZkpt, nrflc, 3);

  printf("\nFolding PXRD reflections back into brillouin zone.\n");  
  /*loop over the number of reflections*/
  for (n=0;n<nrflc;n++) {
    /*find the midpoint of the reciprocal lattice vector defined by the reflection*/
    TTH->BZkpt[n][0] = (double) TTH->rflc_H[n]/2.0;
    TTH->BZkpt[n][1] = (double) TTH->rflc_K[n]/2.0;
    TTH->BZkpt[n][2] = (double) TTH->rflc_L[n]/2.0;
    TTH->hpw[n] = 0;
    TTH->kpw[n] = 0;
    TTH->lpw[n] = 0;
    /*wrap the midpoint into the BZ, while keeping track of how to unwrap it with hklpw*/
    while (TTH->BZkpt[n][0] > 0.5) {
      TTH->hpw[n]++;
      TTH->BZkpt[n][0] = TTH->BZkpt[n][0] - 1;
    }
    while (TTH->BZkpt[n][0] < -0.5) {
      TTH->hpw[n]--;
      TTH->BZkpt[n][0] = TTH->BZkpt[n][0] + 1;
    }
    while (TTH->BZkpt[n][1] > 0.5) {
      TTH->kpw[n]++;
      TTH->BZkpt[n][1] = TTH->BZkpt[n][1] - 1;
    }
    while (TTH->BZkpt[n][1] < -0.5) {
      TTH->kpw[n]--;
      TTH->BZkpt[n][1] = TTH->BZkpt[n][1] + 1;
    }
    while (TTH->BZkpt[n][2] > 0.5) {
      TTH->lpw[n]++;
      TTH->BZkpt[n][2] = TTH->BZkpt[n][2] - 1;
    }
    while (TTH->BZkpt[n][2] < -0.5) {
      TTH->lpw[n]--;
      TTH->BZkpt[n][2] = TTH->BZkpt[n][2] + 1;
    }
  }

} //END of fold_reflection_toBZ 

void symmetry_folded_reflections(TwoTheta * TTH, Symmetry * SYM) 
{
  /*This function finds symmetry related kpts that the reflections
 * can be folded back onto, to minimize number of kpts needed for calculation*/
  int nrflc; /*number of reflections*/
  int nsym; /*number of symmetry elements*/
  int sym; /*symmetry element index*/
  int match; /*integer tracking if a match has been found*/
  int* kpt_match; /*integer tracking if a kpt match has been found*/
  double X_symm, Y_symm, Z_symm; /*symmetry related kpts*/
  int nsym_bzk; /*number of symmetry related kpts*/
  int n, i; /*indices*/
  double* kx_sym; /*symmetry related matches of kpts*/
  double* ky_sym;
  double* kz_sym;
  double pw_x, pw_y, pw_z; /*pws that map the reflections onto the kpts*/
  int k;
  double kx, ky, kz; /*kpt coordinates*/

  /*initialize variables*/
  nrflc = TTH->nrflc;
  nsym = SYM->nsym;

  /*Nullify local arrays*/
  kx_sym = NULL;
  ky_sym = NULL;
  kz_sym = NULL;
  kpt_match = NULL;
  /*Allocate memory for local and structural arrays */
  kpt_match = AllocateMemory_oneD_int(kpt_match, nrflc);
  kx_sym = AllocateMemory_oneD_double(kx_sym, nrflc);
  ky_sym = AllocateMemory_oneD_double(ky_sym, nrflc);
  kz_sym = AllocateMemory_oneD_double(kz_sym, nrflc);
  TTH->kx_sym = AllocateMemory_oneD_double(TTH->kx_sym, nrflc);
  TTH->ky_sym = AllocateMemory_oneD_double(TTH->ky_sym, nrflc);
  TTH->kz_sym = AllocateMemory_oneD_double(TTH->kz_sym, nrflc);
  TTH->BZkpt_sym = AllocateMemory_twoD_double(TTH->BZkpt_sym, nrflc, 3);
  TTH->hpw_sym = AllocateMemory_oneD_int(TTH->hpw_sym, nrflc);
  TTH->kpw_sym = AllocateMemory_oneD_int(TTH->kpw_sym, nrflc);
  TTH->lpw_sym = AllocateMemory_oneD_int(TTH->lpw_sym, nrflc);
  TTH->rflc_H_sym = AllocateMemory_oneD_int(TTH->rflc_H_sym, nrflc);
  TTH->rflc_K_sym = AllocateMemory_oneD_int(TTH->rflc_K_sym, nrflc);
  TTH->rflc_L_sym = AllocateMemory_oneD_int(TTH->rflc_L_sym, nrflc);
  /*end of allocation*/

  printf("\nApplying symmetry to minimize kpts required.\n");

  /*intializing symmetric arrays*/
  kx_sym[0] = TTH->BZkpt[0][0];
  ky_sym[0] = TTH->BZkpt[0][1];
  kz_sym[0] = TTH->BZkpt[0][2];
  nsym_bzk = 1;

  /*Find symmetry related kpts*/
  for (n=0;n<nrflc;n++) {
    match=0;
    for (sym=0;sym<nsym;sym++) {
      /*if there is already a match, break off the search*/
      if (match==1) break;
      /*apply symmetry to one of the stored kpts*/
	  X_symm = (double) SYM->symrel[0][0][sym]*TTH->BZkpt[n][0] + (double) SYM->symrel[1][0][sym]*TTH->BZkpt[n][1] + (double) SYM->symrel[2][0][sym]*TTH->BZkpt[n][2];
	  Y_symm = (double) SYM->symrel[0][1][sym]*TTH->BZkpt[n][0] + (double) SYM->symrel[1][1][sym]*TTH->BZkpt[n][1] + (double) SYM->symrel[2][1][sym]*TTH->BZkpt[n][2];
	  Z_symm = (double) SYM->symrel[0][2][sym]*TTH->BZkpt[n][0] + (double) SYM->symrel[1][2][sym]*TTH->BZkpt[n][1] + (double) SYM->symrel[2][2][sym]*TTH->BZkpt[n][2];
	  
      for (i=0;i<nsym_bzk;i++) {
        /*search through repository of kpts already found, to see if this matches any known kpts*/
		if ((X_symm==kx_sym[i])&&(Y_symm==ky_sym[i])&&(Z_symm==kz_sym[i])) {
          match = 1;
          /*if a match is found, apply the same symmetry to the HKL reflection corresponding to this kpt*/
		  TTH->rflc_H_sym[n] = SYM->symrel[0][0][sym]*TTH->rflc_H[n] + SYM->symrel[1][0][sym]*TTH->rflc_K[n] + SYM->symrel[2][0][sym]*TTH->rflc_L[n];
		  TTH->rflc_K_sym[n] = SYM->symrel[0][1][sym]*TTH->rflc_H[n] + SYM->symrel[1][1][sym]*TTH->rflc_K[n] + SYM->symrel[2][1][sym]*TTH->rflc_L[n];
		  TTH->rflc_L_sym[n] = SYM->symrel[0][2][sym]*TTH->rflc_H[n] + SYM->symrel[1][2][sym]*TTH->rflc_K[n] + SYM->symrel[2][2][sym]*TTH->rflc_L[n];
          kpt_match[n] = i;
          break;
        }
      }

      /*if no match has been found then store this kpt as symmetry-independent one*/
      if (((sym)==(nsym-1))&&(match==0)) {
		kx_sym[nsym_bzk] = TTH->BZkpt[n][0];
		ky_sym[nsym_bzk] = TTH->BZkpt[n][1];
		kz_sym[nsym_bzk] = TTH->BZkpt[n][2];
		TTH->rflc_H_sym[n] = TTH->rflc_H[n];
		TTH->rflc_K_sym[n] = TTH->rflc_K[n];
		TTH->rflc_L_sym[n] = TTH->rflc_L[n];
        kpt_match[n] = nsym_bzk;
		nsym_bzk++;
      }
	}
  }
  TTH->nsym_bzk = nsym_bzk;

  printf("\t Symmetrized kpts and corresponding reflections:\n");
  /*loop over reflections*/
  for (n=0;n<nrflc;n++) {
    /*find midpoints of new reflections*/
    pw_x = (double) TTH->rflc_H_sym[n]/2.0;
    pw_y = (double) TTH->rflc_K_sym[n]/2.0;
    pw_z = (double) TTH->rflc_L_sym[n]/2.0;
    /*find the kpt corresponding to this reflection*/
    k = kpt_match[n];
    kx = kx_sym[k];
    ky = ky_sym[k];
    kz = kz_sym[k];
    /*store these kpts in TTH structure*/
    TTH->BZkpt_sym[n][0] = kx;
    TTH->BZkpt_sym[n][1] = ky;
    TTH->BZkpt_sym[n][2] = kz;
    /*find how these midpoints are unwrapped outside of the 1st BZ*/
    TTH->hpw_sym[n] = pw_x - kx; 
    TTH->kpw_sym[n] = pw_y - ky; 
    TTH->lpw_sym[n] = pw_z - kz; 
    /*print info for each reflection*/
	printf("\t %lf %lf %lf\t%lf %lf %lf %d %d %d\n", pw_x, pw_y, pw_z, TTH->BZkpt_sym[n][0], TTH->BZkpt_sym[n][1], TTH->BZkpt_sym[n][2], TTH->hpw_sym[n], TTH->kpw_sym[n], TTH->lpw_sym[n]);
  }

  /*fill simple symmeterized kpt array*/
  for (i=0;i<nsym_bzk;i++) {
    TTH->kx_sym[i] = kx_sym[i];
    TTH->ky_sym[i] = ky_sym[i];
    TTH->kz_sym[i] = kz_sym[i];
  }

  /*Free local memory*/
  kx_sym = FreeMemory_oneD_double(kx_sym);
  ky_sym = FreeMemory_oneD_double(ky_sym);
  kz_sym = FreeMemory_oneD_double(kz_sym);
  kpt_match = FreeMemory_oneD_int(kpt_match);

} //END of symmetry_folded_reflections

void print_reflections(char filename[200], TwoTheta * TTH, FermiSphere * FS)
{
  /* This function prints the reflection information into a uniform file*/
  FILE * frflc; /*pointer to reflection file*/
  int n; /*index for reflections*/
  int nrflc; /*number of reflections*/
  double intensity; /*intensity of a reflections*/
 
  /*initialize variables*/
  nrflc = TTH->nrflc;
  
  printf("\nPrinting reflection information to: %s.\n", filename);
  /*open the reflection file in write mode*/
  frflc=fopen(filename,"w"); 
  if(frflc==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  /*print header information*/
  fprintf(frflc, "nrflc %d\n", nrflc);
  fprintf(frflc, "FS_angle %lf\n", FS->two_theta);

  /*print powder pattern */
  fprintf(frflc, "\n#__dhkl__2theta__H_K_L__mult__Fhkl_intensity\n");
  for (n=0;n<nrflc;n++) {
    /*find intensity as structure factor times multiplicity*/
    intensity = TTH->F_hkl[n] * TTH->rflc_mult[n];
    fprintf(frflc, "%d\t %lf %lf %d %d %d %d %lf %lf\n", n, TTH->d_hkl[n], TTH->two_theta[n], TTH->rflc_H[n], TTH->rflc_K[n], TTH->rflc_L[n], TTH->rflc_mult[n], TTH->F_hkl[n], intensity);
  }
  
  /*printf symmetrized and folded BZkpts and HKLpw*/
  fprintf(frflc, "\n#__H_K_L___kx_ky_kx____hpw_kpw_lpw\n");
  for  (n=0;n<nrflc;n++) {
    fprintf(frflc, "%d\t%d %d %d\t%lf %lf %lf\t%d %d %d\n", n, TTH->rflc_H_sym[n], TTH->rflc_K_sym[n], TTH->rflc_L_sym[n], TTH->BZkpt_sym[n][0], TTH->BZkpt_sym[n][1], TTH->BZkpt_sym[n][2], TTH->hpw_sym[n], TTH->kpw_sym[n], TTH->lpw_sym[n]);
  }

  /*close reflection file*/
  fclose(frflc);
  
  /*Free allocated variables*/
  TTH->rflc_mult = FreeMemory_oneD_int(TTH->rflc_mult);
  TTH->rflc_H = FreeMemory_oneD_int(TTH->rflc_H);
  TTH->rflc_K = FreeMemory_oneD_int(TTH->rflc_K);
  TTH->rflc_L = FreeMemory_oneD_int(TTH->rflc_L);
  TTH->two_theta = FreeMemory_oneD_double(TTH->two_theta);
  TTH->d_hkl = FreeMemory_oneD_double(TTH->d_hkl);
  TTH->F_hkl = FreeMemory_oneD_double(TTH->F_hkl);
  TTH->hpw = FreeMemory_oneD_int(TTH->hpw);
  TTH->kpw = FreeMemory_oneD_int(TTH->kpw);
  TTH->lpw = FreeMemory_oneD_int(TTH->lpw);
  TTH->BZkpt = FreeMemory_twoD_double(TTH->BZkpt, nrflc);
  TTH->BZkpt_sym = FreeMemory_twoD_double(TTH->BZkpt_sym, nrflc);
  TTH->hpw_sym = FreeMemory_oneD_int(TTH->hpw_sym);
  TTH->kpw_sym = FreeMemory_oneD_int(TTH->kpw_sym);
  TTH->lpw_sym = FreeMemory_oneD_int(TTH->lpw_sym);
  TTH->rflc_H_sym = FreeMemory_oneD_int(TTH->rflc_H_sym);
  TTH->rflc_K_sym = FreeMemory_oneD_int(TTH->rflc_K_sym);
  TTH->rflc_L_sym = FreeMemory_oneD_int(TTH->rflc_L_sym);
} //END of print_reflections

void read_reflections(char filename[200], TwoTheta * TTH)
{
  /*This function reads the reflection file printed from a prep_2theta calculation including the reflections and density associated with each*/
  FILE * frflc; /*pointer to reflection file*/
  int n; 
  int nrflc; /*number of reflections*/
  double twotheta_fs; /*fermi diameter in degrees 2theta*/
  int rflc_no; /*index for reflections*/
  double d_hklv; /*d spacing of reflection*/
  double two_thetav; /*degree 2theta of reflections*/
  int rflc_Hv, rflc_Kv, rflc_Lv; /*HKL indices of reflections*/
  int rflc_multv; /*multiplicity of reflection*/
  double F_hklv; /*structure factor of reflection*/
  double intensity; /*intensity of reflection, squared F_hkl*/
  char strv[50]; /*arbitrary string*/
  double kx, ky, kz; /*reduced kpoint coordinates*/
  int hpwv, kpwv, lpwv; /*hkl indices of planewave*/
  
  printf("\nReading reflection data from: %s.\n", filename);
  /*open the file in read mode*/
  frflc=fopen(filename,"r"); 
  if(frflc==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  
  /*begin reading in reflection data*/

  /*read number of reflections in file*/
  fscanf(frflc, "%s", strv);
  fscanf(frflc, "%d", &nrflc);
  TTH->nrflc = nrflc;

  /*read fermi sphere diameter angle*/
  fscanf(frflc, "%s", strv);
  fscanf(frflc, "%lf", &twotheta_fs);

  /*allocate memroy for reflection information*/
  TTH->rflc_mult = AllocateMemory_oneD_int(TTH->rflc_mult, nrflc);
  TTH->rflc_H = AllocateMemory_oneD_int(TTH->rflc_H, nrflc);
  TTH->rflc_K = AllocateMemory_oneD_int(TTH->rflc_K, nrflc);
  TTH->rflc_L = AllocateMemory_oneD_int(TTH->rflc_L, nrflc);
  TTH->two_theta = AllocateMemory_oneD_double(TTH->two_theta, nrflc);
  TTH->BZkpt = AllocateMemory_twoD_double(TTH->BZkpt, nrflc, 3);
  TTH->hpw = AllocateMemory_oneD_int(TTH->hpw, nrflc);
  TTH->kpw = AllocateMemory_oneD_int(TTH->kpw, nrflc);
  TTH->lpw = AllocateMemory_oneD_int(TTH->lpw, nrflc);
  /*end of allocation*/

  /*read powder pattern */
  fscanf(frflc, "%s", strv);
  for (n=0;n<nrflc;n++) {
    fscanf(frflc, "%d", &rflc_no);
    fscanf(frflc, "%lf", &d_hklv);
    fscanf(frflc, "%lf", &two_thetav);
    fscanf(frflc, "%d %d %d ", &rflc_Hv, &rflc_Kv, &rflc_Lv);
    fscanf(frflc, "%d", &rflc_multv);
    fscanf(frflc, "%lf", &F_hklv);
    fscanf(frflc, "%lf", &intensity);
    /*store necessary variables in TTH struct*/
    TTH->two_theta[n] = two_thetav;
    TTH->rflc_mult[n] = rflc_multv;
  }
  
  /*read symmetrized and folded BZkpts and HKLpw*/
  fscanf(frflc, "%s", strv);
  for  (n=0;n<nrflc;n++) {
    fscanf(frflc, "%d", &rflc_no);
    fscanf(frflc, "%d %d %d", &rflc_Hv, &rflc_Kv, &rflc_Lv);
    fscanf(frflc, "%lf %lf %lf", &kx, &ky, &kz);
    fscanf(frflc, "%d %d %d", &hpwv, &kpwv, &lpwv);
    TTH->rflc_H[n] = rflc_Hv;
    TTH->rflc_K[n] = rflc_Kv;
    TTH->rflc_L[n] = rflc_Lv;
    TTH->BZkpt[n][0] = kx;
    TTH->BZkpt[n][1] = ky;
    TTH->BZkpt[n][2] = kz;
    TTH->hpw[n] = hpwv;
    TTH->kpw[n] = kpwv;
    TTH->lpw[n] = lpwv;
  }

  /*close reflection file*/
  fclose(frflc);

}//END of read_reflections
  
void concatinate_twotheta_potential(EnergyContribution * ECON, EnergyStep * ESTP, TwoTheta *TTH) 
{
  /*This function concatinates the local, nonlocal, and kinetic energy correction terms to find the total potential energy for each refleciton.*/
  int nEstep; /*number of energy bin*/
  int dE; /*energy bin*/
  int nrflc; /*number of reflections*/
  int n; /*reflection index*/
  double reflection_total; /*total potential energy for each reflection*/

  /*initializing variables*/
  nEstep = ESTP->nEstep;
  nrflc = TTH->nrflc;

  /*Allocate Memory*/
  ECON->rflc_total = AllocateMemory_twoD_double(ECON->rflc_total, nEstep, nrflc);
  
  printf("\nConcatinating local and nonlocal potential energy for reflections.\n");
  /*add local and nonlocal grids to find total potential energy*/
  for (n=0;n<nrflc;n++) {
    reflection_total = 0.0;
    for (dE=0;dE<nEstep;dE++) {
      /*add local and nonlocal grids to find total potential energy*/
	  ECON->rflc_total[dE][n] = ECON->rflc_local[dE][n] + ECON->rflc_nonlocal[dE][n] - ECON->rflc_KE_correction[dE][n];
	  reflection_total += ECON->rflc_total[dE][n];
    } /*end dE->nEstep*/
    printf( "\nTotal Potential Energy for rflc %d = %lf \n", n, reflection_total);
  } /*end n->nrflc*/
  
  /*free memory for local, nonlocal, and kecorr grid*/
  ECON->rflc_local = FreeMemory_twoD_double(ECON->rflc_local, nEstep);
  ECON->rflc_nonlocal = FreeMemory_twoD_double(ECON->rflc_nonlocal, nEstep);
  ECON->rflc_KE_correction = FreeMemory_twoD_double(ECON->rflc_KE_correction, nEstep);

} /*end of concatinate_twotheta_potential*/

void search_reflection(char filename[200], MottJonesConditions * MJC, double min_twotheta, double max_twotheta) 
{
  /*This function searches in the reflection file for any 
 * reflections that lie between a min and max 2theta range*/
  FILE * frflc; /*pointer to reflection file*/
  int n; /*reflection index*/
  int nrflc; /*number of reflections*/
  double twotheta_fs; /*fermi sphere diameter in degrees 2theta*/
  int rflc_no; /*reflection index in loop*/
  double d_hklv; /*rec space distance of reflection*/
  double two_thetav; /*two theta of reflection*/
  int rflc_Hv, rflc_Kv, rflc_Lv; /*reflection HKL indices*/
  int* con_H; /*conventional H, K, L indices*/
  int* con_K;
  int* con_L;
  int rflc_multv; /*multiplicity of reflection*/
  double F_hklv; /*structure factor of reflection*/
  double intensity; /*intensity of reflection*/
  char strv[50]; /*arbitrary string*/
  int nfind; /*index for found reflections*/

  printf("\nReading reflection information from: %s\n", filename);
  /*open the filename in read mode*/
  frflc=fopen(filename,"r"); 
  if(frflc==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  /*begin reading header information from rflc file*/
  fscanf(frflc, "%s", strv);
  fscanf(frflc, "%d", &nrflc);
  fscanf(frflc, "%s", strv);
  fscanf(frflc, "%lf", &twotheta_fs);

  /*allocate memory for HKL arrays*/
  con_H = NULL;
  con_K = NULL;
  con_L = NULL;
  con_H = AllocateMemory_oneD_int(con_H, nrflc);
  con_K = AllocateMemory_oneD_int(con_K, nrflc);
  con_L = AllocateMemory_oneD_int(con_L, nrflc);
  
  printf("Find reflections in 2theta range: %lf to %lf\n", min_twotheta, max_twotheta);
  nfind = 0;
  /*read powder pattern */
  for (n=0;n<nrflc;n++) {
    fscanf(frflc, "%d", &rflc_no);
    fscanf(frflc, "%lf", &d_hklv);
    fscanf(frflc, "%lf", &two_thetav);
    fscanf(frflc, "%d %d %d ", &rflc_Hv, &rflc_Kv, &rflc_Lv);
    fscanf(frflc, "%d", &rflc_multv);
    fscanf(frflc, "%lf", &F_hklv);
    fscanf(frflc, "%lf", &intensity);
    
    /*check if reflection is in range; print it if it is*/
    if ((two_thetav >= min_twotheta) && (two_thetav <= max_twotheta)) {
      printf("\t%d\t %lf %lf %d %d %d %d %lf\n", n, d_hklv, two_thetav, rflc_Hv, rflc_Kv, rflc_Lv, rflc_multv, intensity);
      /*convert the HKL indices from primitive to conventional setting*/
      Conventional con = HKL_convert_toP(MJC, rflc_Hv, rflc_Kv, rflc_Lv);
      /*store conventional indices*/
      con_H[nfind] = con.H;      
      con_K[nfind] = con.K;      
      con_L[nfind] = con.L;      
      nfind++;
    }
  }
  
  /*print out conventional indices*/
  printf("\nFound Reflections (conventional):\n");
  for (n=0;n<nfind;n++) {
    printf("\t%d\t %d %d %d\n", n, con_H[n], con_K[n], con_L[n]);
  }

  /*close the reflection file*/
  fclose(frflc);

  /*free local memory*/
  con_H = FreeMemory_oneD_int(con_H);
  con_K = FreeMemory_oneD_int(con_K);
  con_L = FreeMemory_oneD_int(con_L);

}  //END of search_reflection


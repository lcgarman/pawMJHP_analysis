#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "mjhp_HKL_energy.h"
#include "extra_functions.h"

void mjhpHKL_local_energy(NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, VectorIndices *VECT, EnergyStep  *ESTP, EnergyContribution *ECON)
{ 
  /*This function calculate the local potential energy contribution to the mjhp*/
  int nkpt; /*number of kpts*/
  int kptno; /*index for kpts*/
  int nband; /*number of bands*/
  int band; /*index for bands*/
  int npw; /*numper of planewaves*/
  double fermi; /*fermi energy*/
  int ngfftx, ngffty, ngfftz; /*number grid for pws*/
  double occ; /*band occupancy*/
  double wtk; /*kpt weight*/
  double bandE; /*band energy*/
  int pw1, pw2; /*index for planewaves*/
  int h1, k1, l1; /*hkl for periodic part of pw*/
  int h2, k2, l2; /*hkl for periodic part of pw*/
  int delta_h, delta_k, delta_l; /*difference between pw pairs*/
  double pw1_x, pw1_y, pw1_z; /*pw1 reduced coordinates*/
  double pw2_x, pw2_y, pw2_z; /*pw2 reduced coordinates*/
  double ang_pw1x, ang_pw1y, ang_pw1z; /*pw1 coordinates in angstroms*/
  double ang_pw2x, ang_pw2y, ang_pw2z; /*pw2 coordinates in angstroms*/
  double mag_pw1, mag_pw2; /*magnitude of the pws*/
  double mag_diff; /*magnitude of difference bw the pws*/
  int hkl_match; /*tracks if an HKL match has been found*/
  int H_match, K_match, L_match; /*indices of a matching HKL vector*/
  double one_wavecoef_RE, one_wavecoef_IM; /*wavefunction coefficients for pw1*/
  double two_wavecoef_RE, two_wavecoef_IM; /*wavefunction coefficients for pw2*/
  gsl_complex c1, c1_star; /*complex wavefxn coeff for pw1*/
  gsl_complex c2; /*complex wavefxn coeff for pw2*/
  gsl_complex c1star_c2; /*complex conjugate of pw1 coeff times pw2 coeff*/
  gsl_complex Vhkl; /*local potential on reciprocal space grid*/
  gsl_complex complexnums; /*variable storing complex numbers*/
  double noncomplex; /*varibale storing noncomplex numbers*/
  gsl_complex potentialE; /*local energy contribution*/
  double total_local; /*total summed local energy*/

  int j; /*index*/
  int nHKL; /*number of symmetry eq HKL indices*/
  double exponent; /*exponential part of broadening*/
  double sigma; /*extent of broadening parameter*/
  double broad; /*mjhp broadening variable*/
  int nEstep; /*number of band energy bins*/
  int dE; /*band energy bin*/
  double bandE_min; /*minimum band energy*/
  double Emesh; /*density of band energy mesh*/
  double minr, maxr; /*defines the mjhp shell to survey reciprocal space*/
  
  /*initializing variables*/
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi  = WFK->fermi;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;

  nEstep = ESTP->nEstep;
  bandE_min = ESTP->bandE_min;
  minr = VECT->minr;
  maxr = VECT->maxr;
  nHKL = VECT->nHKL;
  sigma = SIGMA;
  Emesh = EMESH;

  /*allocating memory*/
  total_local = 0.0;
  ECON->local = AllocateMemory_oneD_double(ECON->local, nEstep);

  printf( "\nCalculating Local Potential Energy:\n");
  /*start calculating local potential energy contribution for each kpt*/
  for (kptno=0;kptno<nkpt;kptno++){ 
    printf( "kpt %d \t%lf %lf %lf\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno]);
	npw = WFK->npw[kptno];

    /*loop over first planewaves in pair*/
	for(pw1=0;pw1<npw;pw1++) {
      /*store hkl indices for periodic component of pw1*/
	  h1 = WFK->kg[kptno][pw1][0];
	  k1 = WFK->kg[kptno][pw1][1];
	  l1 = WFK->kg[kptno][pw1][2];
      /*find coordinates of pw1 k+G_hkl*/
	  pw1_x = WFK->kpt[0][kptno] + (double) h1;
	  pw1_y = WFK->kpt[1][kptno] + (double) k1;
	  pw1_z = WFK->kpt[2][kptno] + (double) l1;
	  /*find magnitude of pw1 in inverse angstroms*/ 
	  ang_pw1x = pw1_x*UC->ang_ax_star + pw1_y*UC->ang_bx_star + pw1_z*UC->ang_cx_star;
	  ang_pw1y = pw1_x*UC->ang_ay_star + pw1_y*UC->ang_by_star + pw1_z*UC->ang_cy_star;
	  ang_pw1z = pw1_x*UC->ang_az_star + pw1_y*UC->ang_bz_star + pw1_z*UC->ang_cz_star;
	  mag_pw1 = sqrt(ang_pw1x*ang_pw1x+ang_pw1y*ang_pw1y+ang_pw1z*ang_pw1z);
      /*if magnitude of pw1 is outside of mjhp shell go to next iteration*/
      if ((mag_pw1>maxr)||(mag_pw1<minr)) continue;
	  
      /*now loop over second planewave in pair*/
	  for(pw2=0;pw2<npw;pw2++) {
        /*store hkl indices for periodic component of pw2*/
		h2 = WFK->kg[kptno][pw2][0];
		k2 = WFK->kg[kptno][pw2][1];
		l2 = WFK->kg[kptno][pw2][2];
        /*find coordinates of pw2 k+G_hkl*/
		pw2_x = WFK->kpt[0][kptno] + (double) h2;
		pw2_y = WFK->kpt[1][kptno] + (double) k2;
		pw2_z = WFK->kpt[2][kptno] + (double) l2;
	    /*find magnitude of pw2 in inverse angstroms*/ 
		ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
	    mag_pw2 = sqrt(ang_pw2x*ang_pw2x+ang_pw2y*ang_pw2y+ang_pw2z*ang_pw2z);
        /*if magnitude of pw2 is outside of mjhp shell go to next iteration*/
        if ((mag_pw2>maxr)||(mag_pw2<minr)) continue;
		
        /*if pw1 and pw2 both lie in mjhp shell; find difference bw them*/
		delta_h = h1 - h2;
		delta_k = k1 - k2;
		delta_l = l1 - l2;

        /*check symmetry eq HKL to see if delta_hkl matches the indices of interest*/
        hkl_match = 0;
        for (j=0;j<nHKL;j++) {
          if ((delta_h==VECT->H_arr[j])&&(delta_k==VECT->K_arr[j])&&(delta_l==VECT->L_arr[j])) {
            hkl_match = 1;
            H_match = VECT->H_arr[j];
            K_match = VECT->K_arr[j];
            L_match = VECT->L_arr[j];
            break;
		  }
		}
        /*if no match move to next iteration*/
        if (hkl_match == 0) continue;

        /*print out pw matches and coordinates*/
        printf("\tH K L = %d %d %d\n", H_match, K_match, L_match);
        printf("\t\t pw1(%d) = %lf %lf %lf |pw1|=%lf\n", pw1, pw1_x, pw1_y, pw1_z, mag_pw1);
        printf("\t\t pw2(%d) = %lf %lf %lf |pw2|=%lf\n", pw2, pw2_x, pw2_y, pw2_z, mag_pw2);
 
        /*make delta_h positive*/
		if (delta_h < 0) delta_h+=ngfftx;
		if (delta_k < 0) delta_k+=ngffty;
		if (delta_l < 0) delta_l+=ngfftz;

        /*calculate mjhp broadending*/
        mag_diff = mag_pw1 - mag_pw2;
        exponent = -(mag_diff*mag_diff)/sigma;
        broad = exp(exponent);
        printf("\t\tbroadening = %lf\n", broad);

		/*loop over bands to find wavefunction coeff and energy bin dE*/
		for(band=0;band<nband;band++) {
		  /*determine energy bin to store contribution*/
		  bandE = (((WFK->eigen[kptno][band]-fermi)*HATOEV - bandE_min)*Emesh)+0.5;
          dE = floor(bandE);
          /*if energy is above threshold, move on*/
          if ((dE >= nEstep)||(dE < 0)) continue; 
		  
		  /*setting band occupation and kpt weight*/
//		  wtk = 1.0;
//		  occ = WFK->occ[kptno][band];
          wtk = WFK->wtk[kptno];
		  occ = 1.0;
		  
          /*store wavefxn coeffs from pw1 and pw2*/
		  one_wavecoef_RE = WFK->cg[kptno][band][pw1][1];
		  one_wavecoef_IM = WFK->cg[kptno][band][pw1][0];
		  c1 = gsl_complex_rect(one_wavecoef_RE, one_wavecoef_IM);
		  c1_star = gsl_complex_conjugate(c1);
		  two_wavecoef_RE = WFK->cg[kptno][band][pw2][1];
		  two_wavecoef_IM = WFK->cg[kptno][band][pw2][0];
		  c2 = gsl_complex_rect(two_wavecoef_RE, two_wavecoef_IM);
		  
		  /*calculate potential energy for HKL from pw1 and pw2*/
          /* local_E = wtk*occ*(1/cellV)*<c1*|V_HKL|c2> */
		  c1star_c2 = gsl_complex_mul(c1_star, c2);
		  Vhkl = BIN->rec_grid[delta_h][delta_k][delta_l];
		  complexnums = gsl_complex_mul(c1star_c2, Vhkl);
		  noncomplex = occ*wtk*broad/UC->bohr_cellV;
		  potentialE = gsl_complex_mul_real(complexnums, noncomplex);

          /*store local potential energy contribution in structure*/
		  ECON->local[dE] += GSL_REAL(potentialE);
		  total_local += GSL_REAL(potentialE);

		} /*END band->nband loop*/
	  } /*END pw2->npw loop*/
	} /*END pw1->npw loop*/
	printf( "\tkpt %d\t local potential energy = %lf\n", kptno, total_local);
  } /*END: kpt loop*/
  printf( "Total Potential Energy = %lf\n", total_local);

}   //END of mjhpHKL_local_energy function

void mjhpHKL_nonlocal_energy(PawAtomicData * PAW, AtomicVariables * ATM, UnitCell * UC, Wavefunction * WFK, VectorIndices *VECT, EnergyStep  *ESTP, EnergyContribution *ECON)
{
  /*This function calculate the nonlocal potential energy contribution to the mjhp*/
  int nkpt; /*number of kpts*/
  int k; /*kpt index*/
  int nband; /*number of bands*/
  int n; /*band index*/
  int npw; /*numper of planewaves*/
  int pw1, pw2; /*pw indices*/
  int h1, k1, l1; /*hkl for periodic part of pw*/
  int h2, k2, l2; /*hkl for periodic part of pw*/
  double pw1_x, pw1_y, pw1_z; /*reduced coordinates for pw1*/
  double pw2_x, pw2_y, pw2_z; /*reduced coordinates for pw2*/
  double bohr_pw1x, bohr_pw1y, bohr_pw1z; /*coordinates for pw1 in bohr*/
  double bohr_pw2x, bohr_pw2y, bohr_pw2z; /*coordinates for pw1 in bohr*/
  double mag_pw1, mag_pw2; /*magnitude of pws*/
  double g1, g2; /*z of pws in spherical coordinates*/
  double theta1, theta2; /*theta of pws in spherical coordinates*/
  double phi1, phi2; /*phi of pws in spherical coordinates*/
  double costheta1, costheta2; /*cos of the thetas*/
  gsl_complex c1; /*wavefunction coefficient of pw1*/
  gsl_complex c1_star; /*complex conjugate of c1*/
  gsl_complex c2; /*wavefunction coefficient of pw2*/

  double maxr, minr; /*defines the mjhp shell to survey reciprocal space*/
  double ang_pw1x, ang_pw1y, ang_pw1z; /*coordinates for pw1 in ang*/
  double ang_pw2x, ang_pw2y, ang_pw2z; /*coordinates for pw2 in ang*/
  double ang_pw1, ang_pw2; /*magnitude of pws in ang*/
  int delta_h, delta_k, delta_l; /*difference bw pw1 and pw2*/
  int hkl_match; /*tracks if an HKL match has been found*/
  int nHKL; /*number of symmetry equivalent HKL indices*/
  int H_match, K_match, L_match; /*indices of a matching HKL vector*/
  double mag_diff; /*magnitude of difference between pws*/
  double exponent; /*exponent used in mjhp broadening*/
  double sigma; /*extent of broadening parameter*/
  double broad; /*mjhp broading variable*/
  
  int natom; /*number of atoms*/
  int at; /*atom index*/
  int typat; /*typat of atoms*/
  int lnmax; /*number of partial waves*/

  double x, y, z; /*xyz reduced coordinates of atom*/
  double Xcart, Ycart, Zcart; /*xyz cartesian coordinates in bohr*/
  double gg_x, gg_y, gg_z; /*HKL times atom position*/
  double bohr_ggx, bohr_ggy, bohr_ggz; /*HKL times atom position in bohr*/
  double G1G2_R; /*magnitude of HKL times atoms positon in bohr*/
  gsl_complex phaseX; /*phase to account for shift of atom away from origin*/

  int i; /*i->lnm_max keeps track of total l and m values looped over for pw1*/
  int j; /*j->lnm_max keeps track of total l and m values looped over for pw2*/
  int lni, lnj; /*index for partial wave number*/  
  int li, lj; /*angular momentum quantum number for i and j*/

  double fermi; /*fermi energy*/
  double bandE_min; /*minimum band energy*/
  double Emesh; /*energy mesh for band energy*/
  int nEstep; /*number of band energy bins*/
  double occ; /*band occupancy*/
  double wtk; /*kpt weight*/
  double bandE; /*band energy*/
  int dE; /*bin value for band energy*/
  
  double projector_bar2; /*projector bar variable*/
  double projector_bar1; /*projector bar variable*/
  
  int mi, mj; /*magnetic quantum number*/
  int abs_mi, abs_mj; /*absolute value of m*/
  double Plm1, Plm2; /*legendre functions*/
  double Slm1, Slm2; /*real spherical harmonics*/
  double norm; /*normalization factor*/
  double sqrt_two; /*square root of two*/

  double rhoij_RE; /*real component of the rhoij matrix*/
  gsl_complex c1c2; /*wavefunction coefficient*/
  gsl_complex rhoij_IM; /*complex component of rhoij matrix*/
  gsl_complex rhoij; /*rhoij matrix term*/

  double Dij; /*nonlocal matrix element*/
  gsl_complex nonlocalE; /*nonlocal energy*/
  double total_nonlocalE; /*total summed nonlocal energy*/

  double KEij; /*kinetic energy difference matrix element*/
  gsl_complex KE_correction; /*kinetic energy difference*/
  double total_KEcorrection; /*total summed kinetic energy difference*/

  /*initialize variables*/
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi  = WFK->fermi;
  natom = ATM->natom;
  nEstep = ESTP->nEstep;
  bandE_min = ESTP->bandE_min;
  minr = VECT->minr;
  maxr = VECT->maxr;
  nHKL = VECT->nHKL;
  sigma = SIGMA;
  Emesh = EMESH;

  /*precalculate projector_bar variable to save time*/
  int max_pw; /*max number of pws*/
  int max_ln; /*max number of partial waves*/
  double**** projector_bar; /*projector bar variable*/
  projector_bar = NULL;
  /*allocate memory for local projector_bar var*/
  printf("Precalculating projector_bar function...");
  max_pw = 0;
  for (k=0;k<nkpt;k++) {
    if (WFK->npw[k] > max_pw) max_pw = WFK->npw[k];
  }
  max_ln = 0;
  for (at=0;at<natom;at++) {
    typat = ATM->typat[at];
    if (PAW->lnmax[typat] > max_ln) max_ln = PAW->lnmax[typat];
  }
  projector_bar = AllocateMemory_fourD_double(projector_bar, nkpt, max_pw, natom, max_ln);
  projector_bar = calculate_projector_bar(projector_bar, PAW, ATM, WFK, UC);
  printf("Done.\n");
  /*end of projector_bar precalc*/

  /*allocating memory and initializing*/
  total_nonlocalE = 0.0;
  total_KEcorrection = 0.0;
  ECON->nonlocal = AllocateMemory_oneD_double(ECON->nonlocal, nEstep);
  ECON->KE_correction = AllocateMemory_oneD_double(ECON->KE_correction, nEstep);

  /*set constant variables used in loop*/
  norm = (1.0/sqrt(UC->bohr_cellV)) * 4.0*PI;
  sqrt_two = sqrt(2.0);

  printf( "\nCalculating Nonlocal Potential Energy:\n");
  /*start calculating local potential energy contribution for each kpt*/
  for (k=0;k<nkpt;k++) {
    printf( "kpt %d \t%lf %lf %lf\n", k, WFK->kpt[0][k], WFK->kpt[1][k], WFK->kpt[2][k]);
	  
	npw = WFK->npw[k];
    /*loop over first planewaves in pair*/
	for (pw1=0;pw1<npw;pw1++) {
      /*store hkl indices for periodic component of pw1*/
	  h1 = WFK->kg[k][pw1][0];
	  k1 = WFK->kg[k][pw1][1];
	  l1 = WFK->kg[k][pw1][2];
      /*find coordinates of pw1 k+G_hkl*/
	  pw1_x = WFK->kpt[0][k] + (double) h1;
	  pw1_y = WFK->kpt[1][k] + (double) k1;
	  pw1_z = WFK->kpt[2][k] + (double) l1;
	  /*find magnitude of pw1 in inverse bohr*/ 
	  bohr_pw1x = pw1_x*UC->bohr_ax_star + pw1_y*UC->bohr_bx_star + pw1_z*UC->bohr_cx_star;
	  bohr_pw1y = pw1_x*UC->bohr_ay_star + pw1_y*UC->bohr_by_star + pw1_z*UC->bohr_cy_star;
	  bohr_pw1z = pw1_x*UC->bohr_az_star + pw1_y*UC->bohr_bz_star + pw1_z*UC->bohr_cz_star;
	  mag_pw1 = sqrt(bohr_pw1x*bohr_pw1x+bohr_pw1y*bohr_pw1y+bohr_pw1z*bohr_pw1z);
	  /*convert to spherical coordinates*/
	  xyz2sph(bohr_pw1x, bohr_pw1y, bohr_pw1z, &g1, &theta1, &phi1);
	  costheta1 = cos(theta1);
	  
      /*find coordinates in inverse ang for MJHP selection*/
	  ang_pw1x = pw1_x*UC->ang_ax_star + pw1_y*UC->ang_bx_star + pw1_z*UC->ang_cx_star;
	  ang_pw1y = pw1_x*UC->ang_ay_star + pw1_y*UC->ang_by_star + pw1_z*UC->ang_cy_star;
	  ang_pw1z = pw1_x*UC->ang_az_star + pw1_y*UC->ang_bz_star + pw1_z*UC->ang_cz_star;
	  ang_pw1 = sqrt(ang_pw1x*ang_pw1x+ang_pw1y*ang_pw1y+ang_pw1z*ang_pw1z);
      /*if magnitude of pw1 is outside of mjhp shell go to next iteration*/
	  if ((ang_pw1>maxr)||(ang_pw1<minr)) continue;

      /*now loop over second planewave in pair*/
	  for (pw2=0;pw2<npw;pw2++) {
        /*store hkl indices for periodic component of pw2*/
		h2 = WFK->kg[k][pw2][0];
		k2 = WFK->kg[k][pw2][1];
		l2 = WFK->kg[k][pw2][2];
        /*find coordinates of pw2 k+G_hkl*/
		pw2_x = WFK->kpt[0][k] + (double) h2;
		pw2_y = WFK->kpt[1][k] + (double) k2;
		pw2_z = WFK->kpt[2][k] + (double) l2;
	    /*find magnitude of pw2 in inverse bohr*/ 
		bohr_pw2x = pw2_x*UC->bohr_ax_star + pw2_y*UC->bohr_bx_star + pw2_z*UC->bohr_cx_star;
		bohr_pw2y = pw2_x*UC->bohr_ay_star + pw2_y*UC->bohr_by_star + pw2_z*UC->bohr_cy_star;
		bohr_pw2z = pw2_x*UC->bohr_az_star + pw2_y*UC->bohr_bz_star + pw2_z*UC->bohr_cz_star;
		mag_pw2 = sqrt(bohr_pw2x*bohr_pw2x+bohr_pw2y*bohr_pw2y+bohr_pw2z*bohr_pw2z);
		/*convert to spherical coordinates*/
		xyz2sph(bohr_pw2x, bohr_pw2y, bohr_pw2z, &g2, &theta2, &phi2);
		costheta2 = cos(theta2);
		
		/*find coordinates in inverse ang for MJHP selection*/
		ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
		ang_pw2 = sqrt(ang_pw2x*ang_pw2x+ang_pw2y*ang_pw2y+ang_pw2z*ang_pw2z);
        /*if magnitude of pw1 is outside of mjhp shell go to next iteration*/
		if ((ang_pw2>maxr)||(ang_pw2<minr)) continue;
		
        /*if pw1 and pw2 both lie in mjhp shell; find difference bw them*/
		delta_h = h1 - h2;
		delta_k = k1 - k2;
		delta_l = l1 - l2;
		
        /*check symmetry eq HKL to see if delta_hkl matches the indices of interest*/
		hkl_match = 0;
		for (j=0;j<nHKL;j++) {
		  if ((delta_h==VECT->H_arr[j])&&(delta_k==VECT->K_arr[j])&&(delta_l==VECT->L_arr[j])) {
			hkl_match = 1;
			H_match = VECT->H_arr[j];
			K_match = VECT->K_arr[j];
			L_match = VECT->L_arr[j];
			break;
		  }
		}
        /*if no match move to next iteration*/
		if (hkl_match == 0) continue;

		/*print out pw matches and coordinates*/
		printf("\tH K L = %d %d %d\n", H_match, K_match, L_match);
		printf("\t\t pw1(%d) = %lf %lf %lf |pw1|=%lf\n", pw1, pw1_x, pw1_y, pw1_z, ang_pw1);
		printf("\t\t pw2(%d) = %lf %lf %lf |pw2|=%lf\n", pw2, pw2_x, pw2_y, pw2_z, ang_pw2);
		
        /*calculate mjhp broadending*/
		mag_diff = ang_pw1 - ang_pw2;
		exponent = -(mag_diff*mag_diff)/sigma;
		broad = exp(exponent);
		printf("\t\tbroadening = %lf\n", broad);
		
		/*loop over bands to find wavefunction coeff and energy bin dE*/
		for(n=0;n<nband;n++) {
		  /*determine energy bin to store contribution*/
		  bandE = (((WFK->eigen[k][n]-fermi)*HATOEV - bandE_min)*Emesh)+0.5;
		  dE = floor(bandE);
          /*if energy is above threshold, move on*/
		  if ((dE >= nEstep)||(dE < 0)) continue; 
		  
		  /*setting band occupation and kpt weight*/
//		  wtk = 1.0;
//		  occ = WFK->occ[k][band];
          wtk = WFK->wtk[k];
		  occ = 1.0;

          /*store wavefxn coeffs from pw1 and pw2*/
		  c1 = gsl_complex_rect(WFK->cg[k][n][pw1][0], WFK->cg[k][n][pw1][1]);
		  c1_star = gsl_complex_conjugate(c1);
 		  c2 = gsl_complex_rect(WFK->cg[k][n][pw2][0], WFK->cg[k][n][pw2][1]);

          /*begin loop over each of the atoms*/
		  for (at=0;at<natom;at++) {
            typat = ATM->typat[at];
            lnmax = PAW->lnmax[typat]; 

            /*store atomic coordinates*/
		    x = ATM->xred[0][at];
		    y = ATM->xred[1][at];
		    z = ATM->xred[2][at];
		    Xcart = x*UC->bohr_ax+y*UC->bohr_bx+z*UC->bohr_cx; 
		    Ycart = x*UC->bohr_ay+y*UC->bohr_by+z*UC->bohr_cy; 
		    Zcart = x*UC->bohr_az+y*UC->bohr_bz+z*UC->bohr_cz; 

            /*find phase factor to accont for shift from origin*/
            gg_x = (double) h2 - (double) h1;
            gg_y = (double) k2 - (double) k1;
            gg_z = (double) l2 - (double) l1;
		    bohr_ggx = gg_x*UC->bohr_ax_star + gg_y*UC->bohr_bx_star + gg_z*UC->bohr_cx_star;
		    bohr_ggy = gg_x*UC->bohr_ay_star + gg_y*UC->bohr_by_star + gg_z*UC->bohr_cy_star;
		    bohr_ggz = gg_x*UC->bohr_az_star + gg_y*UC->bohr_bz_star + gg_z*UC->bohr_cz_star;
            G1G2_R = bohr_ggx*Xcart + bohr_ggy*Ycart + bohr_ggz*Zcart;
            phaseX = gsl_complex_polar(1.0, G1G2_R);

            /*for each atom loop over the different partial waves*/
            i=0; 
			for (lni=0;lni<lnmax;lni++) {
              li = PAW->l_orbital[typat][lni];

              /*use precalcualted projector_bar*/
              projector_bar1 = projector_bar[k][pw1][at][lni];

			  /*calculate spherical harmonics*/
			  for (mi=(-li);mi<(li+1);mi++) {
				abs_mi = abs(mi);

				/*compute normaljzed legendre polynomjal for spherjcal harmonjcs*/
				Plm1 = gsl_sf_legendre_sphPlm(li, abs_mi, costheta1);
				/*find real spherical harmonics*/
				if (mi==0) Slm1 = 1.0*Plm1;
				else if (mi>0) Slm1 = Plm1 * sqrt_two * cos((double) mi * phi1);
				else if (mi<0) Slm1 = Plm1 * sqrt_two * sin((double) abs_mi * phi1);
				
                /*now loop over j values*/
				j=0; 
				for (lnj=0;lnj<lnmax;lnj++) {
				  lj = PAW->l_orbital[typat][lnj];
				  
                  /*use precalcualted projector_bar*/
                  projector_bar2 = projector_bar[k][pw2][at][lnj];

				  /*now calculate spherical harmonics*/
				  for (mj=(-lj);mj<(lj+1);mj++) {
                    abs_mj = abs(mj);

                    /*compute normaljzed legendre polynomjal for Spherjcal harmonjcs*/
                    Plm2 = gsl_sf_legendre_sphPlm(lj, abs_mj, costheta2);

					/*find real spherical harmonics*/
					if (mj==0) Slm2 = 1.0*Plm2;
					else if (mj>0) Slm2 = Plm2 * sqrt_two * cos((double) mj * phi2);
					else if (mj<0) Slm2 = Plm2 * sqrt_two * sin((double) abs_mj * phi2);
              
					/*compute rhoij[kpt][band] matrix*/
					rhoij_RE = wtk*occ*broad * (norm*norm) * (Slm1*Slm2) * (projector_bar1*projector_bar2); 
					c1c2 = gsl_complex_mul(c1_star, c2);
					rhoij_IM = gsl_complex_mul(c1c2, phaseX);

					rhoij = gsl_complex_mul_real(rhoij_IM, rhoij_RE);
					
                    /*find the nonlocal energy Dij*rhoij*/
					Dij = PAW->Dij_matrix[at][i][j]; 
					nonlocalE = gsl_complex_mul_real(rhoij, Dij);
		            ECON->nonlocal[dE] += GSL_REAL(nonlocalE);
					total_nonlocalE += GSL_REAL(nonlocalE);

                    /*find the kinetic energy correction to the nonlocal energy Kij*rhoij*/
                    if ((li==lj) && (mi==mj)) {
                      /*kronicker delta li,lj and mi,mj*/
					  KEij = PAW->KE_diff_matrix[typat][lni][lnj]; 
					  KE_correction = gsl_complex_mul_real(rhoij, KEij);
                    }
                    else {
                      GSL_REAL(KE_correction) = 0.0;
                      GSL_IMAG(KE_correction) = 0.0;
                    }
		            ECON->KE_correction[dE] += GSL_REAL(KE_correction);
					total_KEcorrection += GSL_REAL(KE_correction);

				    j++;
				  } /*end mj->lj*/
				} /*end lnj->lnmax*/
				i++;
			  } /*end mi->li*/
			} /*end lni->lnmax*/
		  } /*end at->natom*/
		} /*end n->nband*/
	  } /*end pw2->npw*/
	} /*end pw1->npw*/
	printf( "   kpt %d\t nonlocal potential energy = %lf\n", k, total_nonlocalE);
	printf( "   kpt %d\t kinetic energy correction = %lf\n", k, total_KEcorrection);
  } /*end kpt->nkpt*/
  printf( "Total Nonlocal Energy = %lf\n", total_nonlocalE);
  printf( "Total Kinetic Energy Correction = %lf\n", total_KEcorrection);

  /*free projector_bar variable*/
  projector_bar = FreeMemory_fourD_double(projector_bar, nkpt, max_pw, natom);

} //END of mjhpHKL_nonlocal function


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
#include "extra_functions.h"
#include "mjhp_2theta_energy.h"

void mjhp_2theta_local_energy(TwoTheta * TTH, NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, EnergyStep  *ESTP, EnergyContribution *ECON)
{ 
  /*This function calculates the local potential energy contributions for each reflection for the mjhp2theta*/
  int nkpt; /*number of kpts*/
  int kptno; /*kpt index*/
  int nband; /*number of bands*/
  int band; /*band index*/
  int npw; /*numper of planewaves*/
  double fermi; /*fermi energy*/
  int ngfftx, ngffty, ngfftz; /*number grid*/
  double bandE; /*band energy*/
  double occ; /*band occupancy*/
  int pw1, pw2; /*pw pair index*/
  int h1, k1, l1; /*hkl for periodic part of pw1*/
  int h2, k2, l2; /*hkl for periodic part of pw2*/
  int nrflc; /*number of reflections in pxrd*/
  int n; /*reflection index*/
  int delta_h, delta_k, delta_l; /*difference bw planewave pairs*/
  double kx, ky, kz; /*reduced kpt coordinates of unwrapped reflection*/
  int hpw, kpw, lpw; /*unwrapped hkl indices for reflection*/
  int H, K, L; /*HKL indices of the reflection*/
  int rflc_mult; /*multiplicity of reflection*/
  double pw1_x, pw1_y, pw1_z; /*pw1 reduced coordinates*/
  double pw2_x, pw2_y, pw2_z; /*pw2 reduced coordinates*/
  double mag_pw1, mag_pw2; /*magnitude of pws*/
  double one_wavecoef_RE, one_wavecoef_IM; /*wavefunction coefficients*/
  double two_wavecoef_RE, two_wavecoef_IM;
  gsl_complex c1, c1_star; /*wavefunction coefficients concatinated*/
  gsl_complex c2;
  gsl_complex c1star_c2;
  gsl_complex Vhkl; /*local potential on reciprocal space grid*/
  gsl_complex complexnums; /*variable storing complex numbers*/
  double noncomplex; /*varibale storing noncomplex numbers*/
  gsl_complex potentialE; /*local energy contribution*/
  double total_local; /*total summed local energy*/

  int dE; /*band energy bin*/
  int nEstep; /*number of band energy bins*/
  double bandE_min; /*minimum band energy*/
  double Emesh; /*density of band energy mesh*/
  
  /*initializing variables*/
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi  = WFK->fermi;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;

  nEstep = ESTP->nEstep;
  bandE_min = ESTP->bandE_min;
  Emesh = EMESH;
  nrflc = TTH->nrflc;

  /*allocating memeroy and initializing*/
  total_local = 0.0;
  ECON->rflc_local = AllocateMemory_twoD_double(ECON->rflc_local, nEstep, nrflc);

  printf( "\nCalculating Local Potential Energy for Each Reflection:\n");
  /*NOW start Calculating Potential Energy contribution for each REFLECTION*/
  for (n=0;n<nrflc;n++) {
    /*store data for this n reflection*/
    kx = TTH->BZkpt[n][0];  
    ky = TTH->BZkpt[n][1];  
    kz = TTH->BZkpt[n][2];  
    hpw = TTH->hpw[n];
    kpw = TTH->kpw[n];
    lpw = TTH->lpw[n];
    H = TTH->rflc_H[n];
    K = TTH->rflc_K[n];
    L = TTH->rflc_L[n];
    rflc_mult = TTH->rflc_mult[n];
	printf( " reflection# %d\t %d %d %d\n", n, H, K, L);

    /*loop over kpts*/
	for (kptno=0;kptno<nkpt;kptno++){ 
      /*check if coordinates of this kpt match the reflection kpt coordinates*/
	  if ((kx!=WFK->kpt[0][kptno])||(ky!=WFK->kpt[1][kptno])||(kz!=WFK->kpt[2][kptno])) continue;
	  npw = WFK->npw[kptno];

	  /*loop over first planewaves in pair*/
	  for(pw1=0;pw1<npw;pw1++) {
        /*store hkl indices for periodic component of pw1*/
		h1 = WFK->kg[kptno][pw1][0];
		k1 = WFK->kg[kptno][pw1][1];
		l1 = WFK->kg[kptno][pw1][2];
        /*check if pw coordinates match the reflection pw coordinates*/
		if ((h1!=hpw)||(k1!=kpw)||(l1!=lpw)) continue;
        /*find coordinates of pw1 k+G_hkl*/
		pw1_x = WFK->kpt[0][kptno] + (double) h1;
		pw1_y = WFK->kpt[1][kptno] + (double) k1;
		pw1_z = WFK->kpt[2][kptno] + (double) l1;
	    /*find magnitude of pw1 in reduced coordinates*/ 
		mag_pw1 = sqrt(pw1_x*pw1_x+pw1_y*pw1_y+pw1_z*pw1_z);
		
		/*now loop over second planewave in pair*/
		for(pw2=0;pw2<npw;pw2++) {
          /*store hkl indices for periodic component of pw1*/
		  h2 = WFK->kg[kptno][pw2][0];
		  k2 = WFK->kg[kptno][pw2][1];
		  l2 = WFK->kg[kptno][pw2][2];
          /*find coordinates of pw2 k+G_hkl*/
		  pw2_x = WFK->kpt[0][kptno] + (double) h2;
		  pw2_y = WFK->kpt[1][kptno] + (double) k2;
		  pw2_z = WFK->kpt[2][kptno] + (double) l2;
	      /*find magnitude of pw2 in reduced coordinates*/ 
		  mag_pw2 = sqrt(pw2_x*pw2_x+pw2_y*pw2_y+pw2_z*pw2_z);
		
          /*if pw1 data matches reflection, find difference bw pw1 and pw2*/
		  delta_h = h1 - h2;
		  delta_k = k1 - k2;
		  delta_l = l1 - l2;
          /*check if difference is the HKL indices for this reflection*/
		  if ((delta_h!=H)||(delta_k!=K)||(delta_l!=L)) {
			continue;
		  }
          /*check if the magnitudes of pw1 and pw2 are equivalent*/
		  if ((mag_pw1!=mag_pw2)) continue;
          /*print out matches for this reflection*/
		  printf( "\tHKL = %d %d %d\n", delta_h, delta_k, delta_l);
		  printf( "\t\tpw1 = %lf %lf %lf\n", pw1_x, pw1_y, pw1_z);
		  printf( "\t\tpw2 = %lf %lf %lf\n", pw2_x, pw2_y, pw2_z);

          /*make delta_h positive*/
		  if (delta_h < 0) delta_h+=ngfftx;
		  if (delta_k < 0) delta_k+=ngffty;
		  if (delta_l < 0) delta_l+=ngfftz;
		  
		  /*loop over bands to find wavefunction coeff and energy bin dE*/
		  for(band=0;band<nband;band++) {
		    /*determine energy bin to store contribution*/
			bandE = (((WFK->eigen[kptno][band]-fermi)*HATOEV - bandE_min)*Emesh)+0.5;
			dE = floor(bandE);
            /*if energy is above threshold, move on*/
			if ((dE >= nEstep)||(dE < 0)) continue; 
			
		    /*setting band occupation*/
			//occ = WFK->occ[kptno][band];
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
			noncomplex = occ*rflc_mult/UC->bohr_cellV;
			potentialE = gsl_complex_mul_real(complexnums, noncomplex);
			
            /*store potential energy contributions*/
			ECON->rflc_local[dE][n] += GSL_REAL(potentialE);
			total_local += GSL_REAL(potentialE);

		  } /*END band->nband loop*/
		} /*END pw2->npw loop*/
	  } /*END pw1->npw loop*/
  	  printf( "\tkpt %d\t local potential energy = %lf\n", kptno, total_local);
	} /*END: kpt loop*/
  } /*END: rflc loop*/
  printf( "Total Potential Energy = %lf\n", total_local);

} //END of mjhp_2theta_local_energy  

void mjhp_2theta_nonlocal_energy(TwoTheta* TTH, PawAtomicData * PAW, AtomicVariables * ATM, UnitCell * UC, Wavefunction * WFK,EnergyStep  *ESTP, EnergyContribution *ECON)
{
  /*This function calculates the local potential energy contributions for each reflection for the mjhp2theta*/
  int nkpt; /*number of kpts*/
  int k; /*kpt index*/
  int nband; /*number of bands*/
  int n; /*band index*/
  int npw; /*number of planewaves*/
  int pw1, pw2; /*pw index*/
  int h1, k1, l1; /*hkl for periodic part of pw1*/
  int h2, k2, l2; /*hkl for periodic part of pw2*/
  double pw1_x, pw1_y, pw1_z; /*pw1 reduced coordinates*/
  double pw2_x, pw2_y, pw2_z; /*pw2 reduced coordinates*/
  double bohr_pw1x, bohr_pw1y, bohr_pw1z; /*coordinates for pw1 in bohr*/
  double bohr_pw2x, bohr_pw2y, bohr_pw2z; /*coordinates for pw2 in bohr*/
  double mag_pw1, mag_pw2; /*magnitude of pws*/
  double g1, g2; /*z of pws in spherical coordinates*/
  double theta1, theta2; /*theta of pws in spherical coordinates*/
  double phi1, phi2; /*phi of pws in spherical coordinates*/
  double costheta1, costheta2; /*cos of thetas*/
  gsl_complex c1; /*wavefunction coefficients of pw1*/
  gsl_complex c1_star; /*complex conjugate of c1*/
  gsl_complex c2; /*wavefunction coefficients of pw2*/

  int nrflc; /*number of reflections*/
  int nr; /*reflection index*/
  int delta_h, delta_k, delta_l; /*difference bw pw1 and pw2*/
  double kx, ky, kz; /*reduced kpt coordinates of unwrapped reflection*/
  int hpw, kpw, lpw; /*unwrapped hkl indices for reflection*/
  int H, K, L; /*HKL indices of reflection*/
  int rflc_mult; /*multiplicity of reflection*/
  
  int natom; /*number of atoms*/
  int at; /*atom index*/
  int typat; /*type of atom*/
  int lnmax; /*number of partial waves*/

  double x, y, z; /*xyz reduced coordinates of atom*/
  double Xcart, Ycart, Zcart; /*xyz cartesian coordinates in bohr*/
  double gg_x, gg_y, gg_z; /*HKL times atom position*/
  double bohr_ggx, bohr_ggy, bohr_ggz; /*HKL times atom position in bohr*/
  double G1G2_R; /*magnitude of HKL times atom position in bohr*/
  gsl_complex phaseX; /*phase to account for shift of atom away from origin*/

  int i; /*i->lnm_max keeps track of total l and m values looped over for pw1*/
  int j; /*j->lnm_max keeps track of total l and m values looped over for pw2*/
  int lni, lnj; /*index for partial wave number*/  
  int li, lj; /*angular momentum quantum number for i and j*/

  double fermi; /*fermi energy*/
  double bandE_min; /*minimum band energy*/
  double Emesh; /*energy mesh for band energy binning*/
  int nEstep; /*number of band energy bins*/
  double bandE; /*band energy*/
  int dE; /*band energy bin*/
  double occ; /*band occupancy*/
  
  double projector_bar2; /*projector bar variable*/
  double projector_bar1;
  
  int mi, mj; /*magnetic quantum number*/
  int abs_mi, abs_mj; /*absolute value of m*/
  double Plm1, Plm2; /*legendre functions*/
  double Slm1, Slm2; /*real spherical harmonics*/
  double norm; /*normalization factor*/
  double sqrt_two; /*squre root of two*/

  double rhoij_RE; /*real component of the rhoij matrix*/
  gsl_complex c1c2; /*wavefunction coefficients*/
  gsl_complex rhoij_IM; /*complex component of the rhoij matrix*/
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
  nrflc = TTH->nrflc;
  Emesh = EMESH;


  /*precalculate projector_bar variable to save time*/
  int max_pw; /*max number of pws*/
  int max_ln; /*max number of partial waves*/
  double**** projector_bar;
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
  ECON->rflc_nonlocal = AllocateMemory_twoD_double(ECON->rflc_nonlocal, nEstep, nrflc);
  ECON->rflc_KE_correction = AllocateMemory_twoD_double(ECON->rflc_KE_correction, nEstep, nrflc);

  /*set constant variables used in loop*/
  norm = (1.0/sqrt(UC->bohr_cellV)) * 4.0*PI;
  sqrt_two = sqrt(2.0);

  printf( "\nCalculating Nonlocal Potential Energy:\n");
  /*NOW start Calculating Potential Energy contribution for each REFLCTION*/
  for (nr=0;nr<nrflc;nr++) {
    /*store data for this n reflection*/
    kx = TTH->BZkpt[nr][0];  
    ky = TTH->BZkpt[nr][1];  
    kz = TTH->BZkpt[nr][2];  
    hpw = TTH->hpw[nr];
    kpw = TTH->kpw[nr];
    lpw = TTH->lpw[nr];
    H = TTH->rflc_H[nr];
    K = TTH->rflc_K[nr];
    L = TTH->rflc_L[nr];
    rflc_mult = TTH->rflc_mult[nr];
	printf( " reflection# %d\t %d %d %d\n", nr, H, K, L);

    /*loop over kpts*/
	for (k=0;k<nkpt;k++){ 
      /*check if coordinates of this kpt match the reflection kpt coordinates*/
	  if ((kx!=WFK->kpt[0][k])||(ky!=WFK->kpt[1][k])||(kz!=WFK->kpt[2][k])) continue;
	  npw = WFK->npw[k];
	  
	  /*loop over first planewaves in pair*/
	  for (pw1=0;pw1<npw;pw1++) {
        /*store hkl indices for periodic component of pw1*/
		h1 = WFK->kg[k][pw1][0];
		k1 = WFK->kg[k][pw1][1];
		l1 = WFK->kg[k][pw1][2];
        /*check if pw coordinates match the reflection pw coordinates*/
		if ((h1!=hpw)||(k1!=kpw)||(l1!=lpw)) continue;
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
		  
          /*if pw1 data matches reflection, find difference bw pw1 and pw2*/
		  delta_h = h1 - h2;
		  delta_k = k1 - k2;
		  delta_l = l1 - l2;
          /*check if difference is the HKL indices for this reflection*/
		  if ((delta_h!=H)||(delta_k!=K)||(delta_l!=L)) {
			continue;
		  }
          /*check if the magnitudes of pw1 and pw2 are equivalent*/
		  if ((mag_pw1!=mag_pw2)) continue;
          /*print out matches for this reflection*/
		  printf( "\tHKL = %d %d %d\n", delta_h, delta_k, delta_l);
		  printf( "\t\tpw1 = %lf %lf %lf\n", pw1_x, pw1_y, pw1_z);
		  printf( "\t\tpw2 = %lf %lf %lf\n", pw2_x, pw2_y, pw2_z);
		
		  /*loop over bands to find wavefunction coeff and energy bin dE*/
		  for(n=0;n<nband;n++) {
		    /*determine energy bin to store contribution*/
			bandE = (((WFK->eigen[k][n]-fermi)*HATOEV - bandE_min)*Emesh)+0.5;
			dE = floor(bandE);
            /*if energy is above threshold, move on*/
			if ((dE >= nEstep)||(dE < 0)) continue; 
			
			/*setting band occupation*/
			//occ = WFK->occ[kptno][band];
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
			  
              /*find phase factor to account for shift from origin*/
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

					/*calculate spherical harmonics*/
					for (mj=(-lj);mj<(lj+1);mj++) {
					  abs_mj = abs(mj);

					  /*compute normaljzed legendre polynomjal for spherjcal harmonjcs*/
					  Plm2 = gsl_sf_legendre_sphPlm(lj, abs_mj, costheta2);

					  /*find real spherical harmonics*/
					  if (mj==0) Slm2 = 1.0*Plm2;
					  else if (mj>0) Slm2 = Plm2 * sqrt_two * cos((double) mj * phi2);
					  else if (mj<0) Slm2 = Plm2 * sqrt_two * sin((double) abs_mj * phi2);
              
					  /*compute rhoij[kpt][band] matrix*/
					  rhoij_RE = occ*rflc_mult * (norm*norm) * (Slm1*Slm2) * (projector_bar1*projector_bar2); 
					  c1c2 = gsl_complex_mul(c1_star, c2);
					  rhoij_IM = gsl_complex_mul(c1c2, phaseX);

					  rhoij = gsl_complex_mul_real(rhoij_IM, rhoij_RE);
					
					  /*find the total nonlocal energy Dij*rhoij*/
					  Dij = PAW->Dij_matrix[at][i][j]; 
					  nonlocalE = gsl_complex_mul_real(rhoij, Dij);
					  ECON->rflc_nonlocal[dE][nr] += GSL_REAL(nonlocalE);
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
					  ECON->rflc_KE_correction[dE][nr] += GSL_REAL(KE_correction);
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
	} /*kpt->nkpt*/
  } /*nr->nrflc*/
  printf( "Total Nonlocal Energy = %lf\n", total_nonlocalE);
  printf( "Total Kinetic Energy Correction = %lf\n", total_KEcorrection);

  /*free projector_bar variable*/
  projector_bar = FreeMemory_fourD_double(projector_bar, nkpt, max_pw, natom);

} //END of mjhp2theta_nonlocal function


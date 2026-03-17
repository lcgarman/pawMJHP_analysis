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
#include "boxed_band.h"

void Initialize_BoxedBand(BoxedBand * BOX) 
{
  BOX->local = NULL;
  BOX->nonlocal = NULL;
  BOX->kinetic = NULL;
  BOX->KE_correction = NULL;
  BOX->total = NULL;
}

void Box_Allocation(BoxedBand *BOX, NumberGrid * GRD, UnitCell * UC)
{
  int ngfftx, ngffty, ngfftz;
  int maxH, maxK, maxL;
  double max_kx, max_ky, max_kz;
  double max_G;
  int max_Gbin;

  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;

  maxH = ceil(0.5*(double)ngfftx);
  maxK = ceil(0.5*(double)ngffty);
  maxL = ceil(0.5*(double)ngfftz);

  max_kx = maxH*UC->ang_ax_star+maxK*UC->ang_bx_star+maxL*UC->ang_cx_star;
  max_ky = maxH*UC->ang_ay_star+maxK*UC->ang_by_star+maxL*UC->ang_cy_star;
  max_kz = maxH*UC->ang_az_star+maxK*UC->ang_bz_star+maxL*UC->ang_cz_star;

  max_G = sqrt(max_kx*max_kx + max_ky*max_ky + max_kz*max_kz); 
  max_Gbin = ceil(max_G);
  BOX->max_Gbin = max_Gbin;

  BOX->local = AllocateMemory_oneD_double(BOX->local, max_Gbin);
  BOX->nonlocal = AllocateMemory_oneD_double(BOX->nonlocal, max_Gbin);
  BOX->kinetic = AllocateMemory_oneD_double(BOX->kinetic, max_Gbin);
  BOX->KE_correction = AllocateMemory_oneD_double(BOX->KE_correction, max_Gbin);
  BOX->total = AllocateMemory_oneD_double(BOX->total, max_Gbin);
}
 
void print_boxed_band(char filename[200], BoxedBand * BOX)
{
  FILE* fbox;
  int Gbin;
  int nG;
  double total_box;

  nG = BOX->max_Gbin;

  fbox = fopen(filename, "w");
  if(fbox==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  printf("\nPrinting Boxed Band Energy to: %s\n", filename);
  /*print header information*/
  fprintf(fbox, "NoBoxes %d\n", nG);
  fprintf(fbox, "#_local_nonlocal_kinetic_total_KEcorrection_\n");
  fprintf(fbox, "\n");
  /*END of Header Info*/

  /*Print boxed band energy*/
  total_box = 0.0;
  for (Gbin=0;Gbin<nG;Gbin++) {
    BOX->total[Gbin] = BOX->local[Gbin] + BOX->nonlocal[Gbin] + BOX->kinetic[Gbin];
    fprintf(fbox, "%d\t%e\t%e\t%e\t%e\t%e\n", Gbin, BOX->local[Gbin], BOX->nonlocal[Gbin], BOX->kinetic[Gbin], BOX->total[Gbin], BOX->KE_correction[Gbin]);
    total_box += BOX->total[Gbin];
  }
  printf("Total Band Energy = %e\n", total_box);

  /*free energy contribution*/
  BOX->total = FreeMemory_oneD_double(BOX->total);
  BOX->local = FreeMemory_oneD_double(BOX->local);
  BOX->nonlocal = FreeMemory_oneD_double(BOX->nonlocal);
  BOX->KE_correction = FreeMemory_oneD_double(BOX->KE_correction);
  BOX->kinetic = FreeMemory_oneD_double(BOX->kinetic);
} //END of Print_Reflection function 

void boxed_local(NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, BoxedBand * BOX)
{ 
  int nkpt;
  int kptno;
  int nband;
  int band;
  int npw;
  double fermi;
  int ngfftx, ngffty, ngfftz;
  double bandE;
  double occ;
  int pw1, pw2;
  int h1, k1, l1;
  int h2, k2, l2;
  int delta_h, delta_k, delta_l;
  double pw1_x, pw1_y, pw1_z;
  double pw2_x, pw2_y, pw2_z;
  double mag_pw1, mag_pw2;
  double ang_pw1x, ang_pw1y, ang_pw1z;
  double ang_pw2x, ang_pw2y, ang_pw2z;
  double wtk;
  double one_wavecoef_RE, one_wavecoef_IM;
  double two_wavecoef_RE, two_wavecoef_IM;
  gsl_complex c1, c1_star;
  gsl_complex c2;
  gsl_complex c1star_c2;
  gsl_complex Vhkl;
  gsl_complex complexnums;
  double noncomplex;
  gsl_complex potentialE;
  double total_local;
  
  double kx, ky, kz;
  double mag_G;
  int Gbin;
  double Etol;
  Etol = 10.0;

  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi  = WFK->fermi;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  total_local = 0.0;


  printf( "\nCalculating Local Potential Energy:\n");
  /*NOW start Calculating Potential Energy contribution for each kpt*/
  for (kptno=0;kptno<nkpt;kptno++){ 
    printf( "kpt %d \t%lf %lf %lf\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno]);
	npw = WFK->npw[kptno];

    /*loop over first planewaves in pair*/
	for(pw1=0;pw1<npw;pw1++) {
      /*reduced coordinates of planewave*/
	  h1 = WFK->kg[kptno][pw1][0];
	  k1 = WFK->kg[kptno][pw1][1];
	  l1 = WFK->kg[kptno][pw1][2];
	  /*Find magnitude of pw1 in inverse angstroms*/ 
	  pw1_x = WFK->kpt[0][kptno] + (double) h1;
	  pw1_y = WFK->kpt[1][kptno] + (double) k1;
	  pw1_z = WFK->kpt[2][kptno] + (double) l1;
	  ang_pw1x = pw1_x*UC->ang_ax_star + pw1_y*UC->ang_bx_star + pw1_z*UC->ang_cx_star;
	  ang_pw1y = pw1_x*UC->ang_ay_star + pw1_y*UC->ang_by_star + pw1_z*UC->ang_cy_star;
	  ang_pw1z = pw1_x*UC->ang_az_star + pw1_y*UC->ang_bz_star + pw1_z*UC->ang_cz_star;
	  mag_pw1 = sqrt(ang_pw1x*ang_pw1x+ang_pw1y*ang_pw1y+ang_pw1z*ang_pw1z);
	  
      /*now loop over second planewave in pair*/
	  for(pw2=pw1;pw2<npw;pw2++) {
		h2 = WFK->kg[kptno][pw2][0];
		k2 = WFK->kg[kptno][pw2][1];
		l2 = WFK->kg[kptno][pw2][2];
	    /*Find magnitude of pw2*/ 
		pw2_x = WFK->kpt[0][kptno] + (double) h2;
		pw2_y = WFK->kpt[1][kptno] + (double) k2;
		pw2_z = WFK->kpt[2][kptno] + (double) l2;
		ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
	    mag_pw2 = sqrt(ang_pw2x*ang_pw2x+ang_pw2y*ang_pw2y+ang_pw2z*ang_pw2z);
		
        /*if pw1 and pw2 both lie in shell; find the difference bw them*/
		delta_h = h1 - h2;
		delta_k = k1 - k2;
		delta_l = l1 - l2;

//        Gbin = delta_h*delta_h + delta_k*delta_k + delta_l*delta_l;

        /*find the coordinates of each reflection in inverse ang*/
		kx = delta_h*UC->ang_ax_star+delta_k*UC->ang_bx_star+delta_l*UC->ang_cx_star;
		ky = delta_h*UC->ang_ay_star+delta_k*UC->ang_by_star+delta_l*UC->ang_cy_star;
		kz = delta_h*UC->ang_az_star+delta_k*UC->ang_bz_star+delta_l*UC->ang_cz_star;
        /*find magnitude of reflection*/
		mag_G = sqrt(kx*kx + ky*ky + kz*kz);
        if (mag_G == 0.0) Gbin = 0;
        else Gbin = ceil(mag_G);
 
        /*make delta_h positive*/
		if (delta_h < 0) delta_h+=ngfftx;
		if (delta_k < 0) delta_k+=ngffty;
		if (delta_l < 0) delta_l+=ngfftz;

		/*loop over bands to find wavefunction coeff and dE*/
		for(band=0;band<nband;band++) {
		  bandE = (WFK->eigen[kptno][band]-fermi)*HATOEV;
          if ((bandE > Etol) || (bandE <-Etol)) continue;
		  
		  /*Setting Band Occupation and Kpt Weight*/
		  //occ = WFK->occ[kptno][band];
          //wtk = WFK->wtk[kptno];
		  occ = 1.0;
		  wtk = 1.0;
		  
		  /*Calculate potential energy for HKL from pw1 and pw2*/
		  one_wavecoef_RE = WFK->cg[kptno][band][pw1][1];
		  one_wavecoef_IM = WFK->cg[kptno][band][pw1][0];
		  c1 = gsl_complex_rect(one_wavecoef_RE, one_wavecoef_IM);
		  c1_star = gsl_complex_conjugate(c1);
		  two_wavecoef_RE = WFK->cg[kptno][band][pw2][1];
		  two_wavecoef_IM = WFK->cg[kptno][band][pw2][0];
		  c2 = gsl_complex_rect(two_wavecoef_RE, two_wavecoef_IM);
		  
		  c1star_c2 = gsl_complex_mul(c1_star, c2);
		  Vhkl = BIN->rec_grid[delta_h][delta_k][delta_l];
		  complexnums = gsl_complex_mul(c1star_c2, Vhkl);
          if (pw1!=pw2) {
		    noncomplex = 2.0*occ*wtk/UC->bohr_cellV;
          } else  {
		    noncomplex = occ*wtk/UC->bohr_cellV;
          }
		  potentialE = gsl_complex_mul_real(complexnums, noncomplex);
            
		  BOX->local[Gbin] += GSL_REAL(potentialE);

          /*store potential energy contributions*/
		  total_local += GSL_REAL(potentialE);
		} /*END band->nband loop*/
	  } /*END pw2->npw loop*/
	} /*END pw1->npw loop*/
	printf( "   kpt %d\t local potential energy = %lf\n", kptno, total_local);
  } /*END: kpt loop*/
  printf( "Total Potential Energy = %lf\n", total_local);

}   //END of mjhp_hkl_localpot function

void boxed_nonlocal(PawAtomicData * PAW, AtomicVariables * ATM, UnitCell * UC, Wavefunction * WFK, BoxedBand *BOX)
{
  int nkpt;
  int k;
  int nband;
  int n;
  int npw;
  int pw1, pw2;
  int h1, k1, l1;
  int h2, k2, l2;
  double pw1_x, pw1_y, pw1_z;
  double pw2_x, pw2_y, pw2_z;
  double bohr_pw1x, bohr_pw1y, bohr_pw1z;
  double bohr_pw2x, bohr_pw2y, bohr_pw2z;
  double mag_pw1, mag_pw2;
  double g1, g2;
  double theta1, theta2;
  double phi1, phi2;
  double costheta1, costheta2;
  gsl_complex c1;
  gsl_complex c1_star;
  gsl_complex c2;
  int delta_h, delta_k, delta_l;
  
  int natom;
  int at;
  int typat;
  int lnmax;

  double x, y, z;
  double Xcart, Ycart, Zcart;
  double gg_x, gg_y, gg_z;
  double bohr_ggx, bohr_ggy, bohr_ggz;
  double G1G2_R;
  gsl_complex phaseX;

  int i;
  int j;
  int lni, lnj;
  int li, lj;

  double fermi;
  double occ;
  double wtk;
  double bandE;
  
  double projector_bar2;
  double projector_bar1;
  
  int mi, mj;
  int abs_mi, abs_mj;
  double Plm1, Plm2;
  double Slm1, Slm2;
  double norm;
  double sqrt_two;

  double rhoij_RE;
  gsl_complex c1c2;
  gsl_complex rhoij_IM;
  gsl_complex rhoij;

  double Dij;
  gsl_complex nonlocalE;
  double total_nonlocalE;

  double kx, ky, kz;
  double mag_G;
  int Gbin;
  double Etol;
  Etol = 10.0;

  /*initialize variables*/
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi  = WFK->fermi;
  natom = ATM->natom;

  double KEij;
  gsl_complex KE_correction;
  double total_KEcorrection;

  /*precalculate projector_bar variable to save time*/
  int max_pw;
  int max_ln;
  double**** projector_bar;
  projector_bar = NULL;

  /*allocate memory for local projector_bar var*/
  printf("Precalculating projector_bar function");
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
  printf(".\n");
  /*end of projector_bar precalc*/

  /*allocating memory and initializing*/
  total_nonlocalE = 0.0;
  total_KEcorrection = 0.0;
  /*set constant variables used in loop*/
  norm = (1.0/sqrt(UC->bohr_cellV)) * 4.0*PI;
  sqrt_two = sqrt(2.0);
  

  printf( "\nCalculating Nonlocal Potential Energy:\n");
  /*begin loop over kpts*/
  for (k=0;k<nkpt;k++) {
    printf( "kpt %d \t%lf %lf %lf\n", k, WFK->kpt[0][k], WFK->kpt[1][k], WFK->kpt[2][k]);
	  
	npw = WFK->npw[k];
	for (pw1=0;pw1<npw;pw1++) {
	  /*reduced coordinates of planewave*/
	  h1 = WFK->kg[k][pw1][0];
	  k1 = WFK->kg[k][pw1][1];
	  l1 = WFK->kg[k][pw1][2];
	  /*Find magnitude of pw1 in inverse angstroms*/ 
	  pw1_x = WFK->kpt[0][k] + (double) h1;
	  pw1_y = WFK->kpt[1][k] + (double) k1;
	  pw1_z = WFK->kpt[2][k] + (double) l1;
	  /*find K+G vectors in bohr*/
	  bohr_pw1x = pw1_x*UC->bohr_ax_star + pw1_y*UC->bohr_bx_star + pw1_z*UC->bohr_cx_star;
	  bohr_pw1y = pw1_x*UC->bohr_ay_star + pw1_y*UC->bohr_by_star + pw1_z*UC->bohr_cy_star;
	  bohr_pw1z = pw1_x*UC->bohr_az_star + pw1_y*UC->bohr_bz_star + pw1_z*UC->bohr_cz_star;
	  mag_pw1 = sqrt(bohr_pw1x*bohr_pw1x+bohr_pw1y*bohr_pw1y+bohr_pw1z*bohr_pw1z);
	  /*convert |k+G| to spherical coordinates*/
	  xyz2sph(bohr_pw1x, bohr_pw1y, bohr_pw1z, &g1, &theta1, &phi1);
	  costheta1 = cos(theta1);
	  
	  for (pw2=pw1;pw2<npw;pw2++) {
		/*reduced coordinates of planewave*/
		h2 = WFK->kg[k][pw2][0];
		k2 = WFK->kg[k][pw2][1];
		l2 = WFK->kg[k][pw2][2];
		/*Find magnitude of pw1 in inverse angstroms*/ 
		pw2_x = WFK->kpt[0][k] + (double) h2;
		pw2_y = WFK->kpt[1][k] + (double) k2;
		pw2_z = WFK->kpt[2][k] + (double) l2;
		/*find K+G vectors in bohr*/
		bohr_pw2x = pw2_x*UC->bohr_ax_star + pw2_y*UC->bohr_bx_star + pw2_z*UC->bohr_cx_star;
		bohr_pw2y = pw2_x*UC->bohr_ay_star + pw2_y*UC->bohr_by_star + pw2_z*UC->bohr_cy_star;
		bohr_pw2z = pw2_x*UC->bohr_az_star + pw2_y*UC->bohr_bz_star + pw2_z*UC->bohr_cz_star;
		mag_pw2 = sqrt(bohr_pw2x*bohr_pw2x+bohr_pw2y*bohr_pw2y+bohr_pw2z*bohr_pw2z);
		/*convert |k+G| to spherical coordinates*/
		xyz2sph(bohr_pw2x, bohr_pw2y, bohr_pw2z, &g2, &theta2, &phi2);
		costheta2 = cos(theta2);
		
		/*if pw1 and pw2 both lie in shell; find the difference bw them*/
		delta_h = h1 - h2;
		delta_k = k1 - k2;
		delta_l = l1 - l2;
		
//        Gbin = delta_h*delta_h + delta_k*delta_k + delta_l*delta_l;

        /*find the coordinates of each reflection in inverse ang*/
		kx = delta_h*UC->ang_ax_star+delta_k*UC->ang_bx_star+delta_l*UC->ang_cx_star;
		ky = delta_h*UC->ang_ay_star+delta_k*UC->ang_by_star+delta_l*UC->ang_cy_star;
		kz = delta_h*UC->ang_az_star+delta_k*UC->ang_bz_star+delta_l*UC->ang_cz_star;
        /*find magnitude of reflection*/
		mag_G = sqrt(kx*kx + ky*ky + kz*kz);
        if (mag_G == 0) Gbin = 0;
        else Gbin = ceil(mag_G);

		/*continue calculation of nonlocal energy for these planewaves*/
		for (at=0;at<natom;at++) {
		  typat = ATM->typat[at];
		  lnmax = PAW->lnmax[typat]; 

		  /*find phase factor - acconts for shift from origin*/
		  x = ATM->xred[0][at];
		  y = ATM->xred[1][at];
		  z = ATM->xred[2][at];
		  Xcart = x*UC->bohr_ax+y*UC->bohr_bx+z*UC->bohr_cx; 
		  Ycart = x*UC->bohr_ay+y*UC->bohr_by+z*UC->bohr_cy; 
		  Zcart = x*UC->bohr_az+y*UC->bohr_bz+z*UC->bohr_cz; 
		  
		  gg_x = (double) h2 - (double) h1;
		  gg_y = (double) k2 - (double) k1;
		  gg_z = (double) l2 - (double) l1;
		  bohr_ggx = gg_x*UC->bohr_ax_star + gg_y*UC->bohr_bx_star + gg_z*UC->bohr_cx_star;
		  bohr_ggy = gg_x*UC->bohr_ay_star + gg_y*UC->bohr_by_star + gg_z*UC->bohr_cy_star;
		  bohr_ggz = gg_x*UC->bohr_az_star + gg_y*UC->bohr_bz_star + gg_z*UC->bohr_cz_star;
		  G1G2_R = bohr_ggx*Xcart + bohr_ggy*Ycart + bohr_ggz*Zcart;
		  phaseX = gsl_complex_polar(1.0, G1G2_R);
		  
		  i=0; /*i->lnm_max keeps track of total l and m values looped over*/
		  for (lni=0;lni<lnmax;lni++) {
			li = PAW->l_orbital[typat][lni];
			
			/*use precalcualted projector_bar*/
			projector_bar1 = projector_bar[k][pw1][at][lni];
			
			/*now calculate spherical harmonics*/
			for (mi=(-li);mi<(li+1);mi++) {
			  abs_mi = abs(mi);
			
			  /*compute normaljzed legendre polynomjal for Spherjcal harmonjcs*/
			  Plm1 = gsl_sf_legendre_sphPlm(li, abs_mi, costheta1);
			  /*find real spherical harmonics*/
			  if (mi==0) Slm1 = 1.0*Plm1;
			  else if (mi>0) Slm1 = Plm1 * sqrt_two * cos((double) mi * phi1);
			  else if (mi<0) Slm1 = Plm1 * sqrt_two * sin((double) abs_mi * phi1);
				
			  /*now loop over j values*/
			  j=0; /*j->lnm_max keeps track of total l and m values looped over*/
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
				  /*putting band inside pw loops drastically speeds up calc*/
				  for(n=0;n<nband;n++) {
		            bandE = (WFK->eigen[k][n]-fermi)*HATOEV;
                    if ((bandE > Etol) || (bandE <-Etol)) continue;
 
					/*Setting Band Occupation and Kpt Weight*/
					//occ = WFK->occ[k][n];
					//wtk = WFK->wtk[k];
					wtk = 1.0;
					occ = 1.0;
					
					/*find wavefunction coefficients*/
					c1 = gsl_complex_rect(WFK->cg[k][n][pw1][0], WFK->cg[k][n][pw1][1]);
					c1_star = gsl_complex_conjugate(c1);
					c2 = gsl_complex_rect(WFK->cg[k][n][pw2][0], WFK->cg[k][n][pw2][1]);
					
					/*compute rhoij[kpt][band]*/
                    if (pw1!=pw2) {
					  rhoij_RE = 2.0*wtk*occ * (norm*norm) * (Slm1*Slm2) * (projector_bar1*projector_bar2); 
                    } else  {
					  rhoij_RE = wtk*occ * (norm*norm) * (Slm1*Slm2) * (projector_bar1*projector_bar2); 
                    }
					c1c2 = gsl_complex_mul(c1_star, c2);
					rhoij_IM = gsl_complex_mul(c1c2, phaseX);

					rhoij = gsl_complex_mul_real(rhoij_IM, rhoij_RE);
					
                    /*find the total nonlocal energy Dij*rhoij*/
					Dij = PAW->Dij_matrix[at][i][j]; 
					nonlocalE = gsl_complex_mul_real(rhoij, Dij);
		            BOX->nonlocal[Gbin] += GSL_REAL(nonlocalE);
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
		            BOX->KE_correction[Gbin] += GSL_REAL(KE_correction);
					total_KEcorrection += GSL_REAL(KE_correction);

				    }  /*n->nband*/
				    j++;
				  } /*mj=>lj*/
				} /*lnj->lnmax*/
				i++;
			  } /*mi->li*/
			} /*lni->lnmax*/
		  } /*at->natom*/
	  } /*pw2->npw*/
	} /*pw1->npw*/
	printf( "   kpt %d\t nonlocal potential energy = %lf\n", k, total_nonlocalE);
	printf( "   kpt %d\t kinetic energy correction = %lf\n", k, total_KEcorrection);
  } /*k->nkpt*/
  printf( "Total Nonlocal Energy = %lf\n", total_nonlocalE);
  printf( "Total Kinetic Energy Correction = %lf\n", total_KEcorrection);

  /*free projector_bar variable*/
  projector_bar = FreeMemory_fourD_double(projector_bar, nkpt, max_pw, natom);

} //END of mjhpHKL_nonlocal function

void boxed_kinetic(NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BoxedBand *BOX)
{ 
  int nkpt;
  int kptno;
  int nband;
  int band;
  int npw;
  int pw1;
  double occ;
  double wtk;
  double fermi;
  double bandE;
  int h1, k1, l1;
  double pw1_x, pw1_y, pw1_z;
  double bohr_pw1x, bohr_pw1y, bohr_pw1z;
  double pw1x_sq, pw1y_sq, pw1z_sq; 
  double one_wavecoef_RE, one_wavecoef_IM;
  gsl_complex c1;
  double c1_sq;
  double kinetic;
  double total_kinetic;

  int Gbin;
  double Etol;
  Etol = 10.0;
  
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi = WFK->fermi;
  total_kinetic = 0.0;
 
  printf( "\nCalculating Kinetic Energy:\n");
  /*NOW start Calculating Potential Energy contribution for each kpt*/
  for (kptno=0;kptno<nkpt;kptno++){ 
    printf( "kpt %d \t%lf %lf %lf\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno]);
	npw = WFK->npw[kptno];

	/*find |G+k| needed to calc kinetic energy*/
	for(pw1=0;pw1<npw;pw1++) {
	  /*reduced coordinates of planewave - G_hkl*/
	  h1 = WFK->kg[kptno][pw1][0];
	  k1 = WFK->kg[kptno][pw1][1];
	  l1 = WFK->kg[kptno][pw1][2];
	  /*pw1_xyz = k_xyz + G_hkl (rec latt vect Ghkl + wavevector k)*/
	  pw1_x = WFK->kpt[0][kptno] + (double) h1; 
	  pw1_y = WFK->kpt[1][kptno] + (double) k1; 
	  pw1_z = WFK->kpt[2][kptno] + (double) l1; 
	  /*convert to primitive cartesian coordinates*/ 
	  bohr_pw1x = pw1_x*UC->bohr_ax_star + pw1_y*UC->bohr_bx_star + pw1_z*UC->bohr_cx_star;
	  bohr_pw1y = pw1_x*UC->bohr_ay_star + pw1_y*UC->bohr_by_star + pw1_z*UC->bohr_cy_star;
	  bohr_pw1z = pw1_x*UC->bohr_az_star + pw1_y*UC->bohr_bz_star + pw1_z*UC->bohr_cz_star;
	  /*find |G_hkl+k_xyz|^2*/
	  pw1x_sq = bohr_pw1x * bohr_pw1x; 
	  pw1y_sq = bohr_pw1y * bohr_pw1y; 
	  pw1z_sq = bohr_pw1z * bohr_pw1z; 
	  
	  /*loop over bands to find wavefunction coeff and dE*/
	  for(band=0;band<nband;band++) {
		/*Determine Energy range to store contribution*/
		bandE = (WFK->eigen[kptno][band]-fermi)*HATOEV;
        if ((bandE > Etol) || (bandE <-Etol)) continue;
		
		/*Setting Band Occupation and Kpt Weight*/
		//wtk = WFK->wtk[kptno];
		//occ = WFK->occ[kptno][band];
		occ = 1.0;
		wtk = 1.0;
        Gbin = 0;
		
		/*store real and imaginary components of the wavefunction*/
		one_wavecoef_RE = WFK->cg[kptno][band][pw1][1];
		one_wavecoef_IM = WFK->cg[kptno][band][pw1][0];
		c1 = gsl_complex_rect(one_wavecoef_RE, one_wavecoef_IM);
		c1_sq = gsl_complex_abs2(c1);
		/*calculate the kinetic energy*/
		kinetic = 0.5*occ*wtk*c1_sq*(pw1x_sq+pw1y_sq+pw1z_sq);
        BOX->kinetic[Gbin] += kinetic;
		
		/*store kinetic energy contributions*/
		total_kinetic += kinetic;
	  } /*END pw1->npw loop*/
	} /*END band->nband loop*/
    printf( "   kpt %d\t Kinetic Energy = %lf\n", kptno, total_kinetic);
  } /*END: kpt loop*/
  printf("\tTotal Kinetic Energy = %lf\n", total_kinetic);
		
}   //END of kinetic energy function


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "extra_functions.h"

void xyz2sph(double X, double Y, double Z, double * r, double * theta, double * phi)
{
  /*This function converts cartesian coordinates into spherical coordinates. Used in nonlocal energy calc.*/
  /* Z = r cos(theta), X = r sin(theta) cos(phi) Y = r sin(theta)sin(phi) */
  *r = pow((pow(X,2.0)+pow(Y,2.0)+pow(Z,2.0)),0.5);
  if(*r>0.0) {
	*theta=acos(Z/ *r);
	if(fabs(sin(*theta))>0.0) {
	  *phi=acos(X/(*r * sin(*theta)));
	  if(X/(*r * sin(*theta))>1.0) {
		*phi=acos(1.0);
	  }
	  if(X/(*r * sin(*theta))<-1.0) {
		*phi=acos(-1.0);
	  }
	  if(Y<0.0) *phi=-*phi;
	}
	else *phi=0;
  }
  else {
	*phi=0.0;
	*theta=0.0;
  }
}//END of xyz2sph function
        

double**** calculate_projector_bar(double**** projector_bar, PawAtomicData * PAW, AtomicVariables * ATM, Wavefunction * WFK, UnitCell * UC) 
{
  /*This function calculates the integral over the radial grid used in the nonlocal energy*/
  int nkpt; /*number of kpoints*/
  int k; /*kpoint index*/
  int npw; /*number of planewaves*/
  int pw; /*pw index*/
  int h1, k1, l1; /*rec grid indices*/
  double pw_x, pw_y, pw_z; /*pw coordinates*/
  double bohr_pwx, bohr_pwy, bohr_pwz; /*pw coordinates in bohr*/
  double mag_pw; /*magnitude of planewave*/

  int natom; /*number of atoms*/
  int at; /*atom index*/
  int typat; /*type of atom*/
  int lnmax; /*number of partial waves*/
  int ln; /*partial wave index*/
  int l; /*angular momentum quantum number*/
  int mmax; /*number of radial grid points*/
  int r; /*radial grid indice*/
  double r_value; /*radius value of radial grid*/
  double dr; /*derivative of radius value on radial grid*/
  double kG_r; /*planewave times atomic position*/
  double r_projector; /*radius for projector function*/

  /*initialize variables*/
  nkpt = WFK->nkpt;
  natom = ATM->natom;

  /*being calculating the projector bar variable*/
  for (k=0;k<nkpt;k++) {
    npw = WFK->npw[k];
	for (pw=0;pw<npw;pw++) {
      /*store hkl indices for periodic component of pw*/
	  h1 = WFK->kg[k][pw][0];
	  k1 = WFK->kg[k][pw][1];
	  l1 = WFK->kg[k][pw][2];
      /*find coordinates of pw k+G_hkl*/
	  pw_x = WFK->kpt[0][k] + (double) h1;
	  pw_y = WFK->kpt[1][k] + (double) k1;
	  pw_z = WFK->kpt[2][k] + (double) l1;
	  /*find magnitude of pw in inverse bohr*/ 
	  bohr_pwx = pw_x*UC->bohr_ax_star + pw_y*UC->bohr_bx_star + pw_z*UC->bohr_cx_star;
	  bohr_pwy = pw_x*UC->bohr_ay_star + pw_y*UC->bohr_by_star + pw_z*UC->bohr_cy_star;
	  bohr_pwz = pw_x*UC->bohr_az_star + pw_y*UC->bohr_bz_star + pw_z*UC->bohr_cz_star;
	  mag_pw = sqrt(bohr_pwx*bohr_pwx+bohr_pwy*bohr_pwy+bohr_pwz*bohr_pwz);
	  
      /*begin loop over atoms */
      for (at=0;at<natom;at++) {
        /*initialize variables*/
        typat = ATM->typat[at];
        lnmax = PAW->lnmax[typat];
        mmax = PAW->mmax[typat];

        /*loop over the partial waves*/
        for (ln=0;ln<lnmax;ln++) {
          /*store l for each partial wave*/
          l = PAW->l_orbital[typat][ln];

          /*looping over grid points in radial grid*/
		  for (r=0;r<mmax;r++) {
            /*radial value on radial grid*/
			r_value = (PAW->values[typat][r]);
            /*if the radius is less than the paw radius continue with calc*/
			if (r_value <= PAW->paw_rad[typat]) {
			  dr = (PAW->derivatives[typat][r]);
			  kG_r = mag_pw * r_value;
			  r_projector = PAW->projector_function[typat][ln][r];
              /*calculate the projector bar value*/
			  projector_bar[k][pw][at][ln] += r_value*r_value * dr * r_projector * gsl_sf_bessel_jl(l,kG_r);
			} 
		  } /*end of r->mmax*/
		} /*end of ln->lnmax*/
	  } /*end of at->natom*/
	} /*end of pw->npw*/
  } /*end of kpt->nkpt*/

  return projector_bar;

} //END of calculate_projector_bar function


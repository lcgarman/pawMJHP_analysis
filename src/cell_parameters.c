#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "cell_parameters.h"

void Determine_CellParameters(UnitCell* UC, NumberGrid* GRD) 
{  
  /*Calculates real and reciprocal space latt params*/
  printf( "\nCell Parameters:\n");
  printf( "Direct space lattice parameters (bohr): \n");
  printf( "\ta = %lf\t%lf\t%lf\n", UC->bohr_ax, UC->bohr_ay, UC->bohr_az);
  printf( "\tb = %lf\t%lf\t%lf\n", UC->bohr_bx, UC->bohr_by, UC->bohr_bz);
  printf( "\tc = %lf\t%lf\t%lf\n", UC->bohr_cx, UC->bohr_cy, UC->bohr_cz);
  
  /*Determine cell volume in bohr^3 using primitive vectors*/  
  UC->bohr_cellV = UC->bohr_ax*(UC->bohr_by*UC->bohr_cz-UC->bohr_bz*UC->bohr_cy)-\
				   UC->bohr_ay*(UC->bohr_bx*UC->bohr_cz-UC->bohr_bz*UC->bohr_cx)+\
				   UC->bohr_az*(UC->bohr_bx*UC->bohr_cy-UC->bohr_by*UC->bohr_cx);
  printf("\tcell volume (bohr^3) = %lf\n", UC->bohr_cellV);

  /*Determine the voxel volume in bohr^3*/
  UC->voxelV = UC->bohr_cellV/(GRD->ngfftx*GRD->ngffty*GRD->ngfftz);
  printf("\tvoxel volume (bohr^3) = %lf\n", UC->voxelV);

  /*Determine Reciprocal lattice parameters in Bohr*/
  UC->bohr_ax_star = 2*PI*(UC->bohr_by*UC->bohr_cz-UC->bohr_bz*UC->bohr_cy)/UC->bohr_cellV;
  UC->bohr_ay_star = -2*PI*(UC->bohr_bx*UC->bohr_cz-UC->bohr_bz*UC->bohr_cx)/UC->bohr_cellV;
  UC->bohr_az_star = 2*PI*(UC->bohr_bx*UC->bohr_cy-UC->bohr_by*UC->bohr_cx)/UC->bohr_cellV;
  UC->bohr_bx_star = 2*PI*(UC->bohr_cy*UC->bohr_az-UC->bohr_cz*UC->bohr_ay)/UC->bohr_cellV;
  UC->bohr_by_star = -2*PI*(UC->bohr_cx*UC->bohr_az-UC->bohr_cz*UC->bohr_ax)/UC->bohr_cellV;
  UC->bohr_bz_star = 2*PI*(UC->bohr_cx*UC->bohr_ay-UC->bohr_cy*UC->bohr_ax)/UC->bohr_cellV;
  UC->bohr_cx_star = 2*PI*(UC->bohr_ay*UC->bohr_bz-UC->bohr_az*UC->bohr_by)/UC->bohr_cellV;
  UC->bohr_cy_star = -2*PI*(UC->bohr_ax*UC->bohr_bz-UC->bohr_az*UC->bohr_bx)/UC->bohr_cellV;
  UC->bohr_cz_star = 2*PI*(UC->bohr_ax*UC->bohr_by-UC->bohr_ay*UC->bohr_bx)/UC->bohr_cellV;
  //printf("\tReciprocal lattice parameters in (bohr^-1) =\n");
  //printf("\ta* = %lf\t%lf\t%lf\n", UC->bohr_ax_star, UC->bohr_ay_star, UC->bohr_az_star);
  //printf("\tb* = %lf\t%lf\t%lf\n", UC->bohr_bx_star, UC->bohr_by_star, UC->bohr_bz_star);
  //printf("\tc* = %lf\t%lf\t%lf\n", UC->bohr_cx_star, UC->bohr_cy_star, UC->bohr_cz_star);

  /*Convert everything from units bohr to angstorms*/
  UC->ang_ax=UC->bohr_ax*BOHRTOANG;
  UC->ang_ay=UC->bohr_ay*BOHRTOANG;
  UC->ang_az=UC->bohr_az*BOHRTOANG;
  UC->ang_bx=UC->bohr_bx*BOHRTOANG;
  UC->ang_by=UC->bohr_by*BOHRTOANG;
  UC->ang_bz=UC->bohr_bz*BOHRTOANG;
  UC->ang_cx=UC->bohr_cx*BOHRTOANG;
  UC->ang_cy=UC->bohr_cy*BOHRTOANG;
  UC->ang_cz=UC->bohr_cz*BOHRTOANG;
  printf( "Direct space lattice parameters (ang):\n");
  printf( "\ta = %lf\t%lf\t%lf\n", UC->ang_ax, UC->ang_ay, UC->ang_az);
  printf( "\tb = %lf\t%lf\t%lf\n", UC->ang_bx, UC->ang_by, UC->ang_bz);
  printf( "\tc = %lf\t%lf\t%lf\n", UC->ang_cx, UC->ang_cy, UC->ang_cz);

  /*Determine unit cell volume in angstroms*/
  UC->ang_cellV = UC->ang_ax*(UC->ang_by*UC->ang_cz-UC->ang_cy*UC->ang_bz)-\
				  UC->ang_ay*(UC->ang_bx*UC->ang_cz-UC->ang_cx*UC->ang_bz)+\
 				  UC->ang_az*(UC->ang_bx*UC->ang_cy-UC->ang_cx*UC->ang_by);  
  printf("\tcell volume (ang^3) = %lf\n", UC->ang_cellV);
  /*Determine reciprocal lattice parameters in angstroms*/
  UC->ang_ax_star = (2*PI)*(UC->ang_by*UC->ang_cz-UC->ang_cy*UC->ang_bz)/UC->ang_cellV;
  UC->ang_ay_star = -(2*PI)*(UC->ang_bx*UC->ang_cz-UC->ang_cx*UC->ang_bz)/UC->ang_cellV;
  UC->ang_az_star = (2*PI)*(UC->ang_bx*UC->ang_cy-UC->ang_cx*UC->ang_by)/UC->ang_cellV;
  UC->ang_bx_star = (2*PI)*(UC->ang_cy*UC->ang_az-UC->ang_ay*UC->ang_cz)/UC->ang_cellV;
  UC->ang_by_star = -(2*PI)*(UC->ang_cx*UC->ang_az-UC->ang_ax*UC->ang_cz)/UC->ang_cellV;
  UC->ang_bz_star = (2*PI)*(UC->ang_cx*UC->ang_ay-UC->ang_ax*UC->ang_cy)/UC->ang_cellV;
  UC->ang_cx_star = (2*PI)*(UC->ang_ay*UC->ang_bz-UC->ang_by*UC->ang_az)/UC->ang_cellV;
  UC->ang_cy_star = -(2*PI)*(UC->ang_ax*UC->ang_bz-UC->ang_bx*UC->ang_az)/UC->ang_cellV;
  UC->ang_cz_star = (2*PI)*(UC->ang_ax*UC->ang_by-UC->ang_ay*UC->ang_bx)/UC->ang_cellV;
  printf( "Reciprocal lattice parameters (ang^-1) =\n");
  printf( "\ta* = %lf\t%lf\t%lf\n", UC->ang_ax_star, UC->ang_ay_star, UC->ang_az_star);
  printf( "\tb* = %lf\t%lf\t%lf\n", UC->ang_bx_star, UC->ang_by_star, UC->ang_bz_star);
  printf( "\tc* = %lf\t%lf\t%lf\n", UC->ang_cx_star, UC->ang_cy_star, UC->ang_cz_star);
	
}  //END of determine_cellparameters Function

void Calculate_FermiRadius(UnitCell* UC, FermiSphere* FS) 
{
  /*Calculates Fermi radius in reciprocal space distance*/
  double Ghkl_fermi;

  Ghkl_fermi = 2*cbrt(3*PI*PI*FS->vec/UC->ang_cellV);	
  FS->rad_Ef = Ghkl_fermi;
} //END of Calculate_FermiRadius

void Calculate_FermiDegree(UnitCell* UC, FermiSphere* FS) 
{
  /*Calculates Fermi sphere diameter in degrees two theta.*/
  double Ghkl_fermi; /*fermi sphere diameter*/
  double d_hkl; /*reciprocal space distance of FS diameter*/
  double rad_theta; /*conversion of FS diam to degrees radius*/
  double deg_theta; /*conversion to degrees*/
  double two_theta; /*FS diameter in units two theta*/
  double lambda = 0.709319; /*wavelength of Mo Kalpha radiation*/
  
  /*calculate Fermi sphere diameter in degrees 2theta*/
  Ghkl_fermi = 2*cbrt(3*PI*PI*FS->vec/UC->ang_cellV);	
  d_hkl = (2*PI)/Ghkl_fermi;
  rad_theta = asin((1*lambda)/(2*d_hkl));
  deg_theta = rad_theta*(180/PI);
  two_theta = deg_theta*2;
  FS->two_theta = two_theta;
  printf("Fermi sphere diameter = %lf (degrees 2theta)\n", FS->two_theta);
  
} /*end of Calculate_FermiDegree*/


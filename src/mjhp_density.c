#include <fftw3.h>
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

void mjhpHKL_density(NumberGrid *GRD, Wavefunction *WFK, UnitCell * UC, Symmetry *SYM, VectorIndices *VECT, BinaryGrid *BIN)
{
  /*This function calculates the mjhp density associate with 
 * specific HKL vectors and in a defined energy range from the wfks.*/
  int nkpt; /*number of kpts*/
  int kptno; /*kpt index*/
  int nband; /*number of bands*/
  int band; /*band index*/
  int npw; /*number of planewaves*/
  int pw1; /*pw index*/
  int ngfftx, ngffty, ngfftz; /*reciprocal space number grid*/
  int NGX, NGY, NGZ; /*real space number grid*/
  double bandE; /*band energy*/
  double occ; /*band occupation*/
  double fermi; /*fermi energy*/
  double scanE_start, scanE_stop; /*min and max energy to scan*/
  double wtk; /*kpt weight*/
  int h1, k1, l1; /*hkl for periodic part of pw*/
  double pw1_x, pw1_y, pw1_z; /*pw1 reduced coordinates*/
  double pw2_x, pw2_y, pw2_z; /*pw2 reduced coordinates*/
  double ang_pw1x, ang_pw1y, ang_pw1z; /*pw1 coordinates in angstroms*/
  double ang_pw2x, ang_pw2y, ang_pw2z; /*pw2 coordinates in angstroms*/
  double ang_pw1, ang_pw2; /*magnitude of the pws*/
  int MJ_H, MJ_K, MJ_L; /*HKL indices of interest*/
  int nHKL; /*number of symmetry-equivalent HKL indices*/
  int HKL_mult; /*multiplicity of HKL*/
  double minr, maxr; /*defines the mjhp shell to survey reciprocal space*/
  double mag_diff; /*magnitude of diference bw the pws*/
  double exponent; /*exponential part of the broadening*/
  double broad; /*mjhp broadening variable*/
  double sigma; /*extent of broadening parameter*/
  int hpos, kpos, lpos; /*positive hkl indices*/
  gsl_complex c1; /*wavefunction coefficient*/
  double coeff_total; /*summed total of data on grid points*/
  double*** real_grid; /*3D grid in real space the density is stored on*/
  double*** sym_real_grid; /*symmeterized 3D grid in real space*/
  int jx, jy, jz; /*real grid indices*/
  double re_grid; /*real component density for fft*/
  double im_grid; /*imaginary component of density for fft*/
  int i_index; /*input index for fft*/
  int o_index; /*output index for fft*/
  fftw_complex * grid_in; /*fft input grid*/
  fftw_complex * grid_out; /*fft output grid*/
  fftw_plan wfk_den_plan; /*fft plan*/
  int ngfft_size; /*size of nggft grid*/
  int symm_error_flag; /*flag for symmetry grid loops*/
  double zf, yf, xf; /*indices fo symmetry grid loops*/
  double xf2, yf2, zf2; /*indices after symmetr has been applied*/
  int nsym; /*number of symmetry elements*/
  double ix_d, iy_d, iz_d; /*normalized indices after symmetry applied*/
  int stop; /*integer signaling to stop while loop*/
  int ix, iy, iz; 
  int j;

  /*initialize variables*/
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi = WFK->fermi;
  nsym = SYM->nsym;
  minr = VECT->minr;
  maxr = VECT->maxr;
  scanE_start = VECT->scanE_start;
  scanE_stop = VECT->scanE_stop;
  nHKL = VECT->nHKL;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  NGX = ngfftx + 1;
  NGY = ngffty + 1;
  NGZ = ngfftz + 1;
  sigma = SIGMA;
  coeff_total = 0.0;
  real_grid = NULL;
  sym_real_grid = NULL;

  printf("Scanning Energy Range: %lf\t%lf\n", scanE_start, scanE_stop);
  /*Allocating for WavefunctionCoefficients*/
  real_grid = AllocateMemory_threeD_double(real_grid, ngfftx, ngffty, ngfftz);
  sym_real_grid = AllocateMemory_threeD_double(sym_real_grid, ngfftx, ngffty, ngfftz);
  /*Memory Allocation for fftw3*/
  ngfft_size = ngfftx*ngffty*ngfftz;
  grid_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  wfk_den_plan = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, grid_in, grid_out, FFTW_FORWARD, FFTW_MEASURE); 
  if ((grid_in==NULL)||(grid_out==NULL)) {
    printf("ERROR: Memory Allocation Failed for wfk to density fftw\n");
    exit(0);
  }
  /*End of memory allocation*/
  
  printf("\nCalculating HKL electron density from wavefunction.\n");
  /*Begin Calculating Electron density from HKL*/
  for (kptno=0;kptno<nkpt;kptno++) { 
    printf( "kpt %d \t%lf %lf %lf\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno]);
	npw = WFK->npw[kptno];

    /*loop over bands to find wavefunction coeff and energy*/
	for (band=0;band<nband;band++) {
	  /*Scanning energy within minE to maxE range*/
	  bandE = (WFK->eigen[kptno][band]-fermi)*HATOEV;
      if ((bandE > scanE_stop) || (bandE < scanE_start)) continue;
      printf("\tband %d energy = %lf eV\n",  band, bandE);

	  /*zero out grid_in array*/
	  for (h1=0;h1<ngfftx;h1++) {
		for (k1=0;k1<ngffty;k1++) {
		  for (l1=0;l1<ngfftz;l1++) {
			i_index = h1*ngfftz*ngffty+k1*ngfftz+l1;
			grid_in[i_index][REAL] = 0.0;
			grid_in[i_index][IMAG] = 0.0;
		  }
		}
	  }

      /*loop over first pw in pair*/
	  for (pw1=0;pw1<npw;pw1++) {
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
		ang_pw1 = sqrt(ang_pw1x*ang_pw1x+ang_pw1y*ang_pw1y+ang_pw1z*ang_pw1z);
        /*if magnitude of pw1 is outside of mjhp shell go to next iteration*/
		if ((ang_pw1>maxr)||(ang_pw1<minr)) continue;
		
		/*check if this pw could contribute to HKL*/
        HKL_mult = 0;
		for (j=0;j<nHKL;j++) {
		  MJ_H = VECT->H_arr[j];
		  MJ_K = VECT->K_arr[j];
		  MJ_L = VECT->L_arr[j];
		  /*find pw2 = pw1-HKL*/
		  pw2_x = pw1_x - (double) MJ_H;
		  pw2_y = pw1_y - (double) MJ_K;
		  pw2_z = pw1_z - (double) MJ_L;
		  /*convert to ang to compare to MJHP shell*/
		  ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		  ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		  ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
		  ang_pw2 = sqrt(ang_pw2x*ang_pw2x+ang_pw2y*ang_pw2y+ang_pw2z*ang_pw2z);
		  if ((ang_pw2>maxr)||(ang_pw2<minr)) continue;

		  /*if both pw1 and pw2 lie in shell and pw1-pw2=HKL*/
		  HKL_mult++;
		  /*print out pw matches and coordinates*/
		  printf("\tH K L = %d %d %d\n", MJ_H, MJ_K, MJ_L);
		  printf("\t\t pw1(%d) = %lf %lf %lf |pw1|=%lf\n", pw1, pw1_x, pw1_y, pw1_z, ang_pw1);
		  printf("\t\t pw2 = %lf %lf %lf |pw2|=%lf\n", pw2_x, pw2_y, pw2_z, ang_pw2);

		  /*calculate mjhp broadending*/
		  mag_diff = ang_pw1 - ang_pw2;
		  exponent = -(mag_diff*mag_diff)/sigma;
		  broad = exp(exponent);
		  printf("\t\tbroadening = %lf\t", broad);
		  printf("\t\tHKL_mult = %d\n", HKL_mult);
		}
        /*if no possible pw pairing skip this pw*/
		if (HKL_mult == 0) continue;

		/*make hkl grid all (+) numbers*/
		if (h1<0) hpos = h1 + ngfftx;
		else hpos = h1;
		if (k1<0) kpos = k1 + ngffty;
		else kpos = k1;
		if (l1<0) lpos = l1 + ngfftz;
		else lpos = l1;
		i_index = hpos*ngfftz*ngffty+kpos*ngfftz+lpos;

		grid_in[i_index][REAL] = HKL_mult * WFK->cg[kptno][band][pw1][0];
		grid_in[i_index][IMAG] = HKL_mult * -WFK->cg[kptno][band][pw1][1];
	  }
      fftw_execute(wfk_den_plan);

	  /*setting band occupation and kpt weight*/
	  //occ = WFK->occ[kptno][band];
	  occ = 1;
      wtk = WFK->wtk[kptno];
      printf("\t\twtk = %lf\n", wtk);

      /*unwrap fft grid in direct space*/
	  for(jx=0;jx<ngfftx;jx++) {
		for(jy=0;jy<ngffty;jy++) {
		  for(jz=0;jz<ngfftz;jz++) {
			o_index = jx*ngfftz*ngffty+jy*ngfftz+jz;
            re_grid = grid_out[o_index][REAL];
            im_grid = grid_out[o_index][IMAG];
            c1 = gsl_complex_rect(re_grid, im_grid);
            real_grid[jx][jy][jz] += broad*wtk*occ*gsl_complex_abs2(c1)/UC->bohr_cellV;
            coeff_total += broad*wtk*occ*gsl_complex_abs2(c1)/UC->bohr_cellV;

		  }  /*END jz->ngfftz loop*/
		}  /*END jy->ngffty loop*/
	  }  /*END jz->ngfftz loop*/
	} /*END band->nband loop*/
    printf("\tkpt %d\t density = %lf\n", kptno, coeff_total*UC->voxelV);
  } /*END: kpt loop*/
  printf("\nTotal Density = %lf\n", coeff_total*UC->voxelV);

  /*free fftw memory*/
  fftw_destroy_plan(wfk_den_plan);
  fftw_free(grid_in);
  fftw_free(grid_out);

  /*allocate memory for 3d density grid to be stored*/
  BIN->real_grid = AllocateMemory_threeD_double(BIN->real_grid, NGX, NGY, NGZ);

  printf("Symmetrizing the direct space density grid.\n");
  /*symmeterize the density grid in real space*/
  symm_error_flag = 0;
  /*loop over the indices in direct space*/
  for(jx=0;jx<ngfftx;jx++) {
    /*find fractional point on grid*/
	xf = (double)jx/(double)ngfftx;
	for(jy=0;jy<ngffty;jy++) {
	  yf = (double)jy/(double)ngffty;
	  for(jz=0;jz<ngfftz;jz++) {
		zf = (double)jz/(double)ngfftz;
        /*apply symmetry to the grid points*/
	    for(j=0;j<nsym;j++) {
	      xf2 = (double)SYM->symrel[0][0][j]*xf + (double)SYM->symrel[0][1][j]*yf + (double)SYM->symrel[0][2][j]*zf + SYM->tnons[0][j];
	      yf2 = (double)SYM->symrel[1][0][j]*xf + (double)SYM->symrel[1][1][j]*yf + (double)SYM->symrel[1][2][j]*zf + SYM->tnons[1][j];
	      zf2 = (double)SYM->symrel[2][0][j]*xf + (double)SYM->symrel[2][1][j]*yf + (double)SYM->symrel[2][2][j]*zf + SYM->tnons[2][j];
          /*revert coordinates back to non-fractional indices*/
	      ix_d = xf2*(double)ngfftx;
	      iy_d = yf2*(double)ngffty;
	      iz_d = zf2*(double)ngfftz;
          /*convert coordinates to positive values*/
	      stop = 0;
	      while(stop==0) {
	        if(ix_d < 0.0) ix_d += 1.0*ngfftx;
              else stop = 1;
	      }
	      stop = 0;
	      while(stop==0) {
	        if(iy_d < 0.0) iy_d += 1.0*ngffty;
	        else stop = 1;
	      }
	      stop = 0;
	      while(stop==0) {
	        if(iz_d < 0.0) iz_d += 1.0*ngfftz;
	        else stop = 1;
	      }
	      /*round coordinates down to nearest integer*/
	      ix = (int)(floor(ix_d+0.5)); 
	      iy = (int)(floor(iy_d+0.5)); 
	      iz = (int)(floor(iz_d+0.5));
	      /*if the absolute value of the difference between the rounded, 
 * positive coordinate is greater than a threshold throw up an error*/
	      if(fabs((double) ix-ix_d) > 0.001) symm_error_flag = 1;
	      if(fabs((double) iy-iy_d) > 0.001) symm_error_flag = 1;
	      if(fabs((double) iz-iz_d) > 0.001) symm_error_flag = 1;
          /*find remainder of i/ngfft to ensure i is within range (0,ngfft-1)*/
	      ix = ix % ngfftx; 
	      iy = iy % ngffty;
	      iz = iz % ngfftz;
          /*store the symmetrizd direct space grid values*/
	      sym_real_grid[ix][iy][iz] += real_grid[jx][jy][jz]/nsym;
	    }
	  }
    }
  } 
  /*print error if symmetry is not compatible*/
  if(symm_error_flag==1) {
    printf("Warning: ngfft grid spacing is incompatible with space group symmetry.\n");
  }
  /*end of symmetrization*/

  /*store symmeterized values in direct space structure*/
  coeff_total = 0.0;
  for(jx=0;jx<ngfftx;jx++) {
    for(jy=0;jy<ngffty;jy++) {
	  for(jz=0;jz<ngfftz;jz++) {
        /*storing values*/
        BIN->real_grid[jx][jy][jz] = sym_real_grid[jx][jy][jz];
        /*find total summed values*/
        coeff_total += BIN->real_grid[jx][jy][jz];
      }
    }
  }
  printf("Total Normalized Density = %lf\n", coeff_total*UC->voxelV);

  /*rewrap grid so jxyz[max] = jxyz[0]*/
  for(jx=0;jx<NGX;jx++) {
    for(jy=0;jy<NGY;jy++) {
      for(jz=0;jz<NGZ;jz++) {
        if (jx==ngfftx) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][jy][jz];
        if (jy==ngffty) BIN->real_grid[jx][jy][jz] = BIN->real_grid[jx][0][jz];
        if (jz==ngfftz) BIN->real_grid[jx][jy][jz] = BIN->real_grid[jx][jy][0];
        if ((jx==ngfftx)&&(jy==ngffty)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][0][jz];
        if ((jx==ngfftx)&&(jz==ngfftz)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][jy][0];
        if ((jy==ngffty)&&(jz==ngfftz)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[jx][0][0];
        if ((jx==ngfftx)&&(jy==ngffty)&&(jz==ngfftz)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][0][0];
	  }
	}
  }

  /*free grids used only locally*/
  real_grid = FreeMemory_threeD_double(real_grid, ngfftx, ngffty);
  sym_real_grid = FreeMemory_threeD_double(sym_real_grid, ngfftx, ngffty);

} //END of mjhpHKL_density function


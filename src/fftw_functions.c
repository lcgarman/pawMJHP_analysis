#include <fftw3.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "fftw_functions.h"

void FFTon_RealGrid( BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC)
{
  /*This function computes the fast fourier transfrom of data on a grid in direct 
 * space to a grid in reciprocal space */
  int ngfftx, ngffty, ngfftz; /*number grid*/
  int NGX, NGY, NGZ; /*fft number grid*/
  int ngfft_size; /*size of number grid*/
  int i_index; 
  int o_index;
  int jx, jy, jz; /*index for real space grid*/
  int h0, k0, l0; /*index for reciprocal space grid*/
  double coeff_total; /*summed total of data on grid points*/
  double fooo, Fooo; /*direct, reciprocal space origin point*/
  fftw_complex * grid_in; /*fft input grid*/
  fftw_complex * grid_out; /*fft output grid*/
  fftw_plan grid_plan; /*fft plan*/
  double re_rec_grid; /*real component of reciprocal grid*/
  double im_rec_grid; /*imaginary component of reciprocal grid*/

  /*initialize variables*/
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  NGX = GRD->NGX;
  NGY = GRD->NGY;
  NGZ = GRD->NGZ;

  /*Allocating Memory for fftw grids*/
  ngfft_size = ngfftx*ngffty*ngfftz;
  grid_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_plan = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, grid_in, grid_out, FFTW_BACKWARD, FFTW_MEASURE); 
  BIN->rec_grid = AllocateMemory_threeD_complex(BIN->rec_grid, ngfftx, ngffty, ngfftz);
  BIN->cc_rec_grid = AllocateMemory_threeD_complex(BIN->cc_rec_grid, ngfftx, ngffty, ngfftz);

  printf("\nPerforming FFT from direct space to reciprocal space.\n");
  /*fill input fft grid with real space grid data*/
  coeff_total = 0.0;
  for(jx=0;jx<NGX;jx++) {
    for(jy=0;jy<NGY;jy++) {
	  for (jz=0;jz<NGZ;jz++) {
        /*initialize 1D index for grid*/
        i_index = jx*ngfftz*ngffty+jy*ngfftz+jz;
          /*if the index is within ngfft bounds store value*/
          if ((jx<ngfftx)&&(jy<ngffty)&&(jz<ngfftz)) {
		  grid_in[i_index][REAL] = BIN->real_grid[jx][jy][jz];
		  grid_in[i_index][IMAG] = 0.0;
          /*sum the data points on real space grid*/
          coeff_total += grid_in[i_index][REAL];
	    }
	  }
	}
  }

  /*normalize the total summed value*/
  coeff_total = coeff_total*UC->voxelV;
  /*find the normalized origin point in real space*/
  fooo = BIN->real_grid[0][0][0]*UC->bohr_cellV;
  printf("\ttotal values in direct space = %lf\t (f000=%lf)\n", coeff_total, fooo);

  /*execute the fft real->reciprocal*/
  fftw_execute(grid_plan);

  /*unwrap the fft in reciprocal space*/
  coeff_total = 0.0;
  for(h0=0;h0<ngfftx;h0++) {
    for(k0=0;k0<ngffty;k0++) {
      for(l0=0;l0<ngfftz;l0++) {
        /*initialize 1D index for grid*/
        o_index = h0*ngfftz*ngffty+k0*ngfftz+l0;
        /*normalize and store output grid values in reciprocal space*/
	    re_rec_grid = grid_out[o_index][REAL]*UC->voxelV;
	    im_rec_grid = grid_out[o_index][IMAG]*UC->voxelV;
        /*fill reciprocal 3D grid up with values*/
		BIN->rec_grid[h0][k0][l0] = gsl_complex_rect(re_rec_grid, im_rec_grid);
        /*fill reciprocal 3D grid up with complex conjugate values*/
		BIN->cc_rec_grid[h0][k0][l0] = gsl_complex_conjugate(BIN->rec_grid[h0][k0][l0]);
		/*sum the data points on real space grid*/
        coeff_total += re_rec_grid;
	  }
	}
  }
  /*find origin point in reciprocal space*/
  Fooo = GSL_REAL(BIN->rec_grid[0][0][0]);
  printf("\ttotal values in reciprocal space = %lf\t (F000=%lf)\n", coeff_total, Fooo);

  /*Freeing memory used in FFT*/
  fftw_destroy_plan(grid_plan);
  fftw_free(grid_in); fftw_free(grid_out);
  BIN->real_grid = FreeMemory_threeD_double(BIN->real_grid, NGX, NGY);

} //END of FFTon_RealGrid

void FFTon_ReciprocalGrid(BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC, Symmetry* SYM)
{
  /*This function computes the fast fourier transfrom of data on a grid in reciprocal
 * space to a grid in direct space */
  int ngfftx, ngffty, ngfftz; /*number grid*/
  int NGX, NGY, NGZ; /*fft number grid*/
  int ngfft_size; /*size of nggft grid*/
  int h0, k0, l0; /*index for reciprocal space grid*/
  int jx, jy, jz; /*index for direct space grid*/
  int i_index, o_index;
  fftw_complex * grid_in; /*fft input data*/
  fftw_complex * grid_out; /*fft output data*/
  fftw_plan grid_plan; /*fft plan*/
  double re_grid; /*real component of fft grid*/
  double im_grid; /*imaginary component of fft grid*/
  double*** real_grid; /*grid in direct space*/
  double*** sym_real_grid; /*symmetric grid in direct space*/
  double coeff_total; /*total summed values*/
  double Fooo, fooo; /*reciprocal, direct origin point values*/
  int symm_error_flag; /*flag for symmetry mismatches*/
  double zf, yf, xf; /*indices for symmetry grid loops*/
  double xf2, yf2, zf2; /*indices after symmetry has been applied*/
  int nsym; /*number of symmetry elements*/
  double ix_d, iy_d, iz_d; /*normalized indices after symmetry applied*/
  int stop; /*integer signaling to stop while loop*/
  int ix, iy, iz; 
  int j;

  /*initialize variables*/
  nsym = SYM->nsym;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  NGX = ngfftx+1;
  NGY = ngffty+1;
  NGZ = ngfftz+1;
  coeff_total = 0.0;
  real_grid = NULL;
  sym_real_grid = NULL;

  /*allocating memory*/
  real_grid = AllocateMemory_threeD_double(real_grid, ngfftx, ngffty, ngfftz);
  sym_real_grid = AllocateMemory_threeD_double(sym_real_grid, ngfftx, ngffty, ngfftz);
  /*memory allocation for fftw3*/
  ngfft_size = ngfftx*ngffty*ngfftz;
  grid_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_plan = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, grid_in, grid_out, FFTW_FORWARD, FFTW_MEASURE); 
  if ((grid_in==NULL)||(grid_out==NULL)) {
    printf("ERROR: Memory Allocation Failed for wfk to density fftw\n");
    exit(0);
  }
  /*End of memory allocation*/

  printf("\nPerforming FFT from reciprocal space to direct space\n");
  coeff_total = 0.0;
  /*fill fftw 3d grid using values on reciprocal spaced grid*/
  for (h0=0;h0<ngfftx;h0++) {
    for (k0=0;k0<ngffty;k0++) {
      for (l0=0;l0<ngfftz;l0++) {
        /*index for input grid values*/
        i_index = h0*ngfftz*ngffty+k0*ngfftz+l0;
        /*fill fft varibales with reciprocal space grid values*/
		grid_in[i_index][REAL] = GSL_REAL(BIN->rec_grid[h0][k0][l0]);
		grid_in[i_index][IMAG] = GSL_IMAG(BIN->rec_grid[h0][k0][l0]);
        /*sum the values on reciprocal spaced grid*/
		coeff_total += grid_in[i_index][REAL];
      }
    }
  }
  /*find the normalized origin point in reciprocal space*/
  Fooo = GSL_REAL(BIN->rec_grid[0][0][0]);
  printf("\ttotal values in reciprocal space = %lf\t (F000=%lf)\n", coeff_total, Fooo);

  /*execute the fft to bring reciprocal to direct spaced data*/
  fftw_execute(grid_plan);

  coeff_total = 0.0;
  /*unwrap fft in direct space*/
  for (jx=0;jx<ngfftx;jx++) {
    for (jy=0;jy<ngffty;jy++) {
      for (jz=0;jz<ngfftz;jz++) {
        /*index for output grid values*/
        o_index = jx*ngfftz*ngffty+jy*ngfftz+jz;
        /*normalize and store values in direct space*/
        re_grid = grid_out[o_index][REAL]/UC->bohr_cellV;
        im_grid = grid_out[o_index][IMAG]/UC->bohr_cellV;
        /*fill real space grid with coefficients*/
        real_grid[jx][jy][jz] = re_grid*UC->voxelV;
        /*sum the values on reciprocal space grid*/
        coeff_total += real_grid[jx][jy][jz];
      }
    }
  }
  /*normalize and store origin point in direct space*/
  fooo = real_grid[0][0][0]*ngfft_size;
  printf("\ttotal values in direct space = %lf (f000 = %lf)\n", coeff_total, fooo);

  /*free memory needed for fft calculation*/
  fftw_destroy_plan(grid_plan);
  fftw_free(grid_in);
  fftw_free(grid_out);
  BIN->rec_grid = FreeMemory_threeD_complex(BIN->rec_grid, ngfftx, ngffty);
  
  /*allocate memory for real space potential grid*/
  BIN->real_grid = AllocateMemory_threeD_double(BIN->real_grid, NGX, NGY, NGZ);

  printf("\tSymmetrizing the direct space grid.\n");
  /*symmeterize the potential grid in real space*/
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
	      /* round coordinates down to nearest integer*/
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
          /*store the symmetrized direct space grid values*/
	      sym_real_grid[ix][iy][iz] += real_grid[jx][jy][jz]/nsym;
	    }
	  }
    }
  }
  /*print error if symmetry is not compatible*/
  if(symm_error_flag==1) {
    printf("Warning: ngfft grid spacing is incompatible with space group symmetry.\n");
  }
  /*end of symmeterization*/

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
  /*find normalized origin point values*/
  fooo = BIN->real_grid[0][0][0]*ngfft_size;
  printf("\ttotal symmetrized value =  %lf (f000=%lf)\n", coeff_total, fooo);

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

} //END of FFTon_ReciprocalGrid function



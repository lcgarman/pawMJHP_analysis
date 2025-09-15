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
#include "cell_parameters.h"
#include "fftw_functions.h"
#include "paw_parameters.h"
#include "extra_functions.h"
#include "local_energy.h"
#include "nonlocal_energy.h"
#include "kinetic_energy.h"
#include "kinetic_correction.h"
#include "mjhp_functions.h"
#include "hkl_functions.h"
#include "energy_functions.h"
#include "rwa_functions.h"
#include "read_binary_v10.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char MJOUTfilename [200];
  char MJOUTfilename_N [200];
  char ABOfilename [200];
  char ABOUTfilename [200];
  char WFKfilename [200];
  char POTfilename [200];
  char H_append [20];
  char K_append [20];
  char L_append [20];
  int H_conventional, K_conventional, L_conventional;
  int n;
  int ndts;

  /*declaring structres*/
  FileCabinet fcab;
  UnitCell ucell;
  NumberGrid grid; 
  Symmetry symm;
  Wavefunction wave;
  BinaryGrid bin;
  AtomicVariables atom;
  PawAtomicData paw;
  EnergyStep estp;
  EnergyContribution econ;
  VectorIndices vect;
  FermiSphere fsph;
  MottJonesConditions mjc;

  /*copy argument string to inupt filename*/
  if (argc > 1) {
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: paw_mjhpHKL <*mjin>  ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_PawAtomicData(&paw);
  Initialize_EnergyContribution(&econ);
  Initialize_VectorIndices(&vect);
  Initialize_MottJonesConditions(&mjc);

  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);
  ndts = mjc.ndts;
  
  /*Store filenames */
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJOUTfilename, fcab.MJOUTfilename);  
  printf("input: %s\n", ABOfilename);
  printf("output: %s\n", MJOUTfilename);
  printf("\nNumber of HKL: %d\n", ndts);
  for (n=0;n<ndts;n++) {
    printf("\tDTSET %d\t HKL = %d %d %d\n", n+1, mjc.jzH[n], mjc.jzK[n], mjc.jzL[n]);
  }
  printf("\nBegin Reading Files.\n");

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  printf("\nwavefunction: %s\n", WFKfilename);
  read_binary_v10(WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);

  /*find energy bounds for mjhpHKL calc*/
  find_energy_bounds(&wave, &estp);

  /*calculate potential energy*/
  strcpy(POTfilename, ABOfilename);
  strcat(POTfilename, "_o_POT");
  printf("\npotential: %s\n", POTfilename);
  read_binary_v10(POTfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 
  /*fft pot grid to reciprocal space*/
  FFTon_RealGrid(&bin, &grid, &ucell);
  /*free memory for 3d grid that will not be used*/
  bin.cc_rec_grid = FreeMemory_threeD_complex(bin.cc_rec_grid, grid.ngfftx, grid.ngffty);

  /*read in paw psp values from abinit output file*/
  strcpy(ABOUTfilename, ABOfilename);
  strcat(ABOUTfilename, ".out");
  printf("\nout_file: %s\n", ABOUTfilename);
  read_about_pseudopotential(ABOUTfilename, &paw, &atom);    
  read_atomic_data(&paw, &atom);

  /*Begin calculating potential for datasets*/
  for (n=0;n<ndts;n++) {
    /*store HKL indices as conventional*/
    mjc.nH = mjc.jzH[n];
    mjc.nK = mjc.jzK[n];
    mjc.nL = mjc.jzL[n];
    H_conventional = mjc.nH;
    K_conventional = mjc.nK;
    L_conventional = mjc.nL;
    printf("\nDTSET %d\t HKL = %d %d %d\n", n+1, mjc.nH, mjc.nK, mjc.nL);

    /*transform HKL to primitive centering*/
    transform_HKL(&mjc);
    /*store primitive HKL in VectorIndices structure*/
    vect.H = mjc.nH;
    vect.K = mjc.nK;
    vect.L = mjc.nL;
    printf("\nHKL Conventional: %d %d %d\n", H_conventional, K_conventional, L_conventional);
    printf("HKL Primitive: %d %d %d\n", vect.H, vect.K, vect.L);

    /*find nonsymmorphic symmetry related HKL vectors*/
    find_symmetric_hkl(&vect, &symm, &grid); 

    /*find the shell in rec space to include pot energy contributions from*/
    find_MJregion(&vect, &ucell);

    /*Calculate the Local potential energy contribution*/
    mjhpHKL_local_energy(&grid, &wave, &ucell, &bin, &vect, &estp, &econ);
  
    /*calculate the nonlocal potential energy contribution*/
    mjhpHKL_nonlocal_energy(&paw, &atom, &ucell, &wave, &vect, &estp, &econ);
    
    /*combine nonlocal and local potential energy*/
    concatinate_HKL_potential(&econ, &estp); 
  
    /*Integrate the total potential energy up to Ef*/
    integrate_HKL_potential(&econ, &estp, &atom);

	/*print the total potential energy contributions*/
	sprintf(H_append, "%d", H_conventional);
	sprintf(K_append, "%d", K_conventional);
	sprintf(L_append, "%d", L_conventional);
	strcpy(MJOUTfilename_N, MJOUTfilename);
	strcat(MJOUTfilename_N,"_");
	strcat(MJOUTfilename_N, H_append);
	strcat(MJOUTfilename_N, K_append);
	strcat(MJOUTfilename_N, L_append);
	strcat(MJOUTfilename_N,".mjout");
	print_mjhpHKL_energy(MJOUTfilename_N, &vect, &estp, &ucell, &econ); 
	
	/*Free allocated memory for next iteration*/ 
	vect.H_arr = FreeMemory_oneD_int(vect.H_arr);
	vect.K_arr = FreeMemory_oneD_int(vect.K_arr);
	vect.L_arr = FreeMemory_oneD_int(vect.L_arr);
	vect.H_posarr = FreeMemory_oneD_int(vect.H_posarr);
	vect.K_posarr = FreeMemory_oneD_int(vect.K_posarr);
	vect.L_posarr = FreeMemory_oneD_int(vect.L_posarr);
  }
  
  /*Free allocated memory*/
  printf( "\nFree Allocated Variables.\n");
  FreeMemory_PAWvariable(&paw, &atom); 
  FreeMemory_Wavefunctions(&wave);
  bin.rec_grid = FreeMemory_threeD_complex(bin.rec_grid, grid.ngfftx, grid.ngffty);
  bin.cc_rec_grid= FreeMemory_threeD_complex(bin.cc_rec_grid, grid.ngfftx, grid.ngffty);
  wave.npw = FreeMemory_oneD_int(wave.npw);
  symm.symrel = FreeMemory_threeD_int(symm.symrel, 3, 3);
  atom.typat = FreeMemory_oneD_int(atom.typat);
  wave.kpt = FreeMemory_twoD_double(wave.kpt, 3);
  symm.tnons = FreeMemory_twoD_double(symm.tnons, 3);
  wave.wtk = FreeMemory_oneD_double(wave.wtk);
  symm.mult = FreeMemory_oneD_int(symm.mult);
  atom.xred = FreeMemory_twoD_double(atom.xred, 3);
  wave.eigen = FreeMemory_twoD_double(wave.eigen, wave.nkpt);
  mjc.jzH = FreeMemory_oneD_int(mjc.jzH);
  mjc.jzK = FreeMemory_oneD_int(mjc.jzK);
  mjc.jzL = FreeMemory_oneD_int(mjc.jzL);
  mjc.scanE_min = FreeMemory_oneD_double(mjc.scanE_min);
  mjc.scanE_mid = FreeMemory_oneD_double(mjc.scanE_mid);
  mjc.scanE_max = FreeMemory_oneD_double(mjc.scanE_max);

  printf( "\n\nEND OF FILE\n");

  return(0);
}

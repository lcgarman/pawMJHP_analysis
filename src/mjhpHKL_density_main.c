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
#include "read_binary_v10.h"
#include "rwa_functions.h"
#include "xsf_functions.h"
#include "cell_parameters.h"
#include "fftw_functions.h"
#include "mjhp_functions.h"
#include "hkl_functions.h"
#include "mjhp_density.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char DENXSF [200];
  char DENXSF_N [200];
  char DENXSF_N1 [200];
  char DENXSF_N2 [200];
  char ABOfilename [200];
  char WFKfilename [200];
  char H_append [20];
  char K_append [20];
  char L_append [20];
  char scanEmin_append [20];
  char scanEmid_append [20];
  char scanEmax_append [20];
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
  EnergyContribution econ;
  VectorIndices vect;
  FermiSphere fsph;
  MottJonesConditions mjc;

  /*copy argument string to input filename*/
  if (argc > 1) {
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: paw_mjhpHKL_density <*mjin>  ");
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
  strcpy(DENXSF, fcab.MJOUTfilename);  
  printf("input: %s\n", ABOfilename);
  printf("output xsf: %s\n", DENXSF);
  printf("\n Number of HKL: %d\n", ndts);
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
    printf("\nHKL Primitive: %d %d %d\n", vect.H, vect.K, vect.L);

    /*find nonsymmorphic symmetry related HKL vectors*/
    find_symmetric_hkl(&vect, &symm, &grid); 

    /*find the shell in rec space to include density contributions from*/
    find_MJregion(&vect, &ucell);

    /*rename the density.xsf files with the energy range and HKL indices*/
	sprintf(H_append, "%d", H_conventional);
	sprintf(K_append, "%d", K_conventional);
	sprintf(L_append, "%d", L_conventional);
	sprintf(scanEmin_append, "%0.3f", mjc.scanE_min[n]);
	sprintf(scanEmid_append, "%0.3f", mjc.scanE_mid[n]);
	sprintf(scanEmax_append, "%0.3f", mjc.scanE_max[n]);
	strcpy(DENXSF_N, DENXSF);
	strcat(DENXSF_N,"_");
	strcat(DENXSF_N, H_append);
	strcat(DENXSF_N, K_append);
	strcat(DENXSF_N, L_append);
	strcat(DENXSF_N,"_");
	strcpy(DENXSF_N1, DENXSF_N);
	strcpy(DENXSF_N2, DENXSF_N);
	strcat(DENXSF_N1, scanEmin_append);
	strcat(DENXSF_N1, "_");
	strcat(DENXSF_N1, scanEmid_append);
	strcat(DENXSF_N1, ".xsf");
	strcat(DENXSF_N2, scanEmid_append);
	strcat(DENXSF_N2, "_");
	strcat(DENXSF_N2, scanEmax_append);
	strcat(DENXSF_N2, ".xsf");

    /*find the density associated with the HKL vectors LOWER SPLIT*/
    vect.scanE_start = mjc.scanE_min[n];
    vect.scanE_stop = mjc.scanE_mid[n];
	mjhpHKL_density(&grid, &wave, &ucell, &symm, &vect, &bin);
	printf( "\nPrinting MJHP HKL Densities to: %s\n", DENXSF_N1);
    print_XSF(DENXSF_N1, &ucell, &grid, &bin, &atom);

    /*find the density associated with the HKL vectors HIGHER SPLIT*/
    vect.scanE_start = mjc.scanE_mid[n];
    vect.scanE_stop = mjc.scanE_max[n];
	mjhpHKL_density(&grid, &wave, &ucell, &symm, &vect, &bin);
	printf( "\nPrinting MJHP HKL Densities to: %s\n", DENXSF_N2);
    print_XSF(DENXSF_N2, &ucell, &grid, &bin, &atom);

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

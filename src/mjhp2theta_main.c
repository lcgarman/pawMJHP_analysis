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
#include "cell_parameters.h"
#include "rwa_functions.h"
#include "fftw_functions.h"
#include "paw_parameters.h"
#include "extra_functions.h"
#include "energy_functions.h"
#include "twotheta_functions.h"
#include "mjhp_2theta_energy.h"
#include "hkl_functions.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char MJOUTfilename [200];
  char RFLCfilename [200];
  char ABOfilename [200];
  char ABOUTfilename [200];
  char POTfilename [200];
  char DENfilename [200];
  char WFKfilename [200];
  double fermi;

  /*declaring structres*/
  FileCabinet fcab;
  UnitCell ucell;
  NumberGrid grid; 
  EnergyStep estp;
  FermiSphere fsph; 
  Symmetry symm;
  Wavefunction wave;
  BinaryGrid bin;
  AtomicVariables atom;
  PawAtomicData paw;
  MottJonesConditions mjc;
  TwoTheta tth;
  EnergyContribution econ;

  /*copy argument string to inupt filename*/
  if (argc > 1) {
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: paw_mjhp2theta <*.mjin> ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_PawAtomicData(&paw);
  Initialize_MottJonesConditions(&mjc);
  Initialize_TwoTheta(&tth);
  Initialize_EnergyContribution(&econ);
  
  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJOUTfilename, fcab.MJOUTfilename);
  printf( "input: %s\n", ABOfilename);
  printf( "output: %s\n", MJOUTfilename);
  
  /*read reflection file*/
  strcpy(RFLCfilename, MJOUTfilename);
  strcat(RFLCfilename, ".rflc");
  printf( "\nrflc: %s\n", RFLCfilename);
  read_reflections(RFLCfilename, &tth); 

  /*Read COMP_i_DEN for Fermi energy file*/
  strcpy(DENfilename, ABOfilename);
  strcat(DENfilename, "_i_DEN");
  printf("\ninput density: %s\n", DENfilename);
  read_binary_v10(DENfilename, 2, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 
  fermi = wave.fermi;

  /*Read potential file*/
  strcpy(POTfilename, ABOfilename);
  strcat(POTfilename, "_o_POT");
  printf("\npotential: %s\n", POTfilename);
  read_binary_v10(POTfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom, &paw);

  /*fine unit cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);
  
  /*Perform FFT on Real space potential grid*/
  FFTon_RealGrid(&bin, &grid, &ucell);

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  printf("\nwavefunction: %s\n", WFKfilename);
  read_binary_v10(WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 
  printf("\n");

  /*Find band energy range to scan*/
  wave.fermi = fermi; /*restores fermi energy saved before from*/
  find_energy_bounds(&wave, &estp); 

  /*calculate the local potential energy contribution*/
  mjhp_2theta_local_energy(&tth, &grid, &wave, &ucell, &bin, &estp, &econ);

  /*read in paw psp values from abinit output file*/
  strcpy(ABOUTfilename, ABOfilename);
  strcat(ABOUTfilename, ".out");
  printf("\nout_file: %s\n", ABOUTfilename);
  read_about_pseudopotential(ABOUTfilename, &paw, &atom);    
  read_atomic_data(&paw, &atom);

  /*calculate the nonlocal potential energy contribution*/
  mjhp_2theta_nonlocal_energy(&tth, &paw, &atom, &ucell, &wave, &estp, &econ);

  /*combine nonlocal and local potential energy*/
  concatinate_twotheta_potential(&econ, &estp, &tth); 
  
  /*calculate the Fermi sphere angle*/
  Calculate_FermiDegree(&ucell, &fsph);
 
  /*print total potential energy contributions*/ 
  strcat(MJOUTfilename, "_2theta");
  strcat(MJOUTfilename, ".mjout");
  print_mjhp2theta_energy(MJOUTfilename, &tth, &estp, &ucell, &econ, &fsph);

  /*Free allocated memory not needed anymore*/
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

  tth.rflc_mult = FreeMemory_oneD_int(tth.rflc_mult);
  tth.two_theta = FreeMemory_oneD_double(tth.two_theta);
  tth.BZkpt = FreeMemory_twoD_double(tth.BZkpt, tth.nrflc);
  tth.hpw = FreeMemory_oneD_int(tth.hpw);
  tth.kpw = FreeMemory_oneD_int(tth.kpw);
  tth.lpw = FreeMemory_oneD_int(tth.lpw);
  tth.rflc_H = FreeMemory_oneD_int(tth.rflc_H);
  tth.rflc_K = FreeMemory_oneD_int(tth.rflc_K);
  tth.rflc_L = FreeMemory_oneD_int(tth.rflc_L);
  mjc.jzH = FreeMemory_oneD_int(mjc.jzH);
  mjc.jzK = FreeMemory_oneD_int(mjc.jzK);
  mjc.jzL = FreeMemory_oneD_int(mjc.jzL);
  mjc.scanE_min = FreeMemory_oneD_double(mjc.scanE_min);
  mjc.scanE_mid = FreeMemory_oneD_double(mjc.scanE_mid);
  mjc.scanE_max = FreeMemory_oneD_double(mjc.scanE_max);

  printf( "\n\nEND OF FILE\n");

  return(0);
}

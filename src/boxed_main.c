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
#include "boxed_band.h"
#include "rwa_functions.h"
#include "read_binary_v10.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char BOXfilename [200];
  char ABOfilename [200];
  char ABOUTfilename [200];
  char WFKfilename [200];
  char POTfilename [200];

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
  BoxedBand box;

  if (argc > 1) {
      //copy file string to inupt files
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: box_the_band <*mjin>  ");
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
  Initialize_BoxedBand(&box);

  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);
  
  /*Store filenames */
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(BOXfilename, fcab.MJOUTfilename);  
  printf("in: %s\n", ABOfilename);
  printf("out: %s\n", BOXfilename);
  printf("Begin Reading Files...\n");

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  printf("wavefunction: %s\n", WFKfilename);
  read_binary_v10(WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);

  /*calculate potential energy*/
  strcpy(POTfilename, ABOfilename);
  strcat(POTfilename, "_o_POT");
  printf("potential: %s\n", POTfilename);
  read_binary_v10(POTfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 
  /*fft pot grid to reciprocal space*/
  FFTon_RealGrid(&bin, &grid, &ucell);
  bin.cc_rec_grid = FreeMemory_threeD_complex(bin.cc_rec_grid, grid.ngfftx, grid.ngffty);

  /*read in paw psp values from abinit output file*/
  strcpy(ABOUTfilename, ABOfilename);
  strcat(ABOUTfilename, ".out");
  printf("out_file: %s\n", ABOUTfilename);
  read_about_pseudopotential(ABOUTfilename, &paw, &atom);    
  read_atomic_data(&paw, &atom);

  Initialize_BoxedBand(&box); 
  Box_Allocation(&box, &grid, &ucell);

  boxed_local(&grid, &wave, &ucell, &bin, &box);

  boxed_nonlocal(&paw, &atom, &ucell, &wave, &box);

  boxed_kinetic(&grid, &wave, &ucell, &box);

  strcat(BOXfilename, "_box");
  print_boxed_band(BOXfilename, &box);

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

  printf( "\n\nEND OF FILE\n");

  return(0);
}

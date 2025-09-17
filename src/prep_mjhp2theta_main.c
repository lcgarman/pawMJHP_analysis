#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include "globals.h"
#include "structures.h"
#include "read_binary_v10.h"
#include "allocate_memory.h"
#include "rwa_functions.h"
#include "cell_parameters.h"
#include "fftw_functions.h"
#include "recgrid_functions.h"
#include "twotheta_functions.h"
#include "hkl_functions.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char RFLCfilename [200];
  char ABOfilename [200];
  char ABINfilename [200];
  char MABINfilename [200];
  char DENfilename [200];

  /*declaring structres*/
  FileCabinet fcab;
  UnitCell ucell;
  NumberGrid grid; 
  FermiSphere fsph;
  Symmetry symm;
  Wavefunction wave;
  BinaryGrid bin;
  AtomicVariables atom;
  MottJonesConditions mjc;
  TwoTheta tth;
  PawAtomicData paw;

  /*copy argument string to input filename*/
  if (argc > 1) {
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: prep_paw2theta <*.mjin> ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_MottJonesConditions(&mjc);
  Initialize_TwoTheta(&tth);
  Initialize_PawAtomicData(&paw);
  
  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(RFLCfilename, fcab.MJOUTfilename);
  printf("input: %s\n", ABOfilename);
  printf("output: %s\n", RFLCfilename);

  /*Read input density file*/
  strcpy(DENfilename, ABOfilename);
  strcat(DENfilename, "_i_DEN");
  printf("\ninput density: %s\n", DENfilename);
  read_binary_v10(DENfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom, &paw); 

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);

  /*calculate the Fermi sphere angle*/
  Calculate_FermiDegree(&ucell, &fsph);

  /*define the reciprocal space grid*/
  find_HKLgrid_max(&grid); 
  prepare_HKLgrid(&grid);
  
  /*Perform FFT on Real space density grid*/
  FFTon_RealGrid(&bin, &grid, &ucell);

  /*calculate the simulated powder pattern*/
  calculate_powder_pattern(&tth, &bin, &grid, &ucell, &symm); 
  /*fold the reflections back into the 1st BZ*/
  fold_reflections_toBZ(&tth);
  /*apply symmetry to minimize the number of kpts needed to analyze*/
  symmetry_folded_reflections(&tth, &symm);

  /*print the reflection information*/
  strcat(RFLCfilename, ".rflc");
  print_reflections(RFLCfilename, &tth, &fsph);

  /*read abinit compound.in file and modify input variables*/
  strcpy(ABINfilename, ABOfilename);
  strcat(ABINfilename, ".in");
  strcpy(MABINfilename, ABINfilename);
  strcat(MABINfilename, "-modified");
  modify_abinitin(ABINfilename, MABINfilename, &tth);
  /*rename modified abinit in file to original in filename*/
  remove(ABINfilename);
  rename(MABINfilename, ABINfilename);

  /*Free allocated memory not needed anymore*/
  printf( "\nFree Allocated Variables.\n");
  FreeMemory_PAWvariable(&paw, &atom); 
  bin.rec_grid = FreeMemory_threeD_complex(bin.rec_grid, grid.ngfftx, grid.ngffty);
  bin.cc_rec_grid= FreeMemory_threeD_complex(bin.cc_rec_grid, grid.ngfftx, grid.ngffty);
  wave.npw = FreeMemory_oneD_int(wave.npw);
  symm.symrel = FreeMemory_threeD_int(symm.symrel, 3, 3);
  atom.typat = FreeMemory_oneD_int(atom.typat);
  wave.kpt = FreeMemory_twoD_double(wave.kpt, 3);
  symm.tnons = FreeMemory_twoD_double(symm.tnons, 3);
  symm.mult = FreeMemory_oneD_int(symm.mult);
  atom.xred = FreeMemory_twoD_double(atom.xred, 3);

  printf("\n\nEND OF FILE\n");

  return(0);
}

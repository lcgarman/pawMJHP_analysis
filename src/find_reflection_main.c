#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "globals.h"
#include "rwa_functions.h"
#include "allocate_memory.h"
#include "twotheta_functions.h"
#include "hkl_functions.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char RFLCfilename [200];
  char MJINfilename [200];
  double min_twotheta; 
  double max_twotheta;

  /*declaring structures*/
  FileCabinet fcab;
  FermiSphere fsph;
  MottJonesConditions mjc;
  TwoTheta tth;

  /*if two arguments given copy arg to input filename and manually define theta range*/
  if (argc == 2) {
    strcpy(MJINfilename, argv[1]);
    min_twotheta = 0.0;
    max_twotheta =  60.0;
  }
  /*if four arguments given copy arg to input filename and store the min and max theta*/
  else if (argc == 4) {
    strcpy(MJINfilename, argv[1]);
    min_twotheta = atof(argv[2]);
    max_twotheta = atof(argv[3]);
  }
  else {
    printf("USAGE: find_reflections <*.mjin> ");
	exit(0);
  }

  /*initialzie all pointers in struct to NULL*/
  Initialize_MottJonesConditions(&mjc);
  Initialize_TwoTheta(&tth);

  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);

  /*store filenames*/
  strcpy(RFLCfilename, fcab.MJOUTfilename);
  strcat(RFLCfilename, ".rflc");

  /*search reflections for inside min to max 2theta*/
  search_reflection(RFLCfilename, &mjc, min_twotheta, max_twotheta);

  /*Free allocated variables*/
  mjc.jzH = FreeMemory_oneD_int(mjc.jzH);
  mjc.jzK = FreeMemory_oneD_int(mjc.jzK);
  mjc.jzL = FreeMemory_oneD_int(mjc.jzL);
  mjc.scanE_min = FreeMemory_oneD_double(mjc.scanE_min);
  mjc.scanE_mid = FreeMemory_oneD_double(mjc.scanE_mid);
  mjc.scanE_max = FreeMemory_oneD_double(mjc.scanE_max);

  return(0);

}

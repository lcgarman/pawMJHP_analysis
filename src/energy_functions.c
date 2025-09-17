#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"

void find_energy_bounds(Wavefunction * WFK, EnergyStep * ESTP) 
{
  /*This function bins the band energy, and finds the max and min band energy values*/
  int k; /*kpt index*/
  int n; /*band index*/
  int nkpt; /*number of kpts*/
  int nband; /*number of bands*/
  double fermi; /*fermi energy*/
  double bandE_min_Ha; /*minimum band energy*/
  double bandE_max_Ha; /*maximum band energy*/
  int bandE_min_integer; /*min band energy rounded*/
  int bandE_max_integer; /*max band energy rounded*/
  double zero; /*zero of band energy*/

  /*initializing variables*/ 
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi = WFK->fermi;

  /*find the minimum and maximum band energy*/
  printf( "\nFinding energy range to scan.\n");
  bandE_min_Ha = 0.0;
  bandE_max_Ha = 0.0;
  for (k=0;k<nkpt;k++) {
    for (n=0;n<nband;n++) {
      if (WFK->eigen[k][n] < bandE_min_Ha) bandE_min_Ha = WFK->eigen[k][n];
      if (WFK->eigen[k][n] > bandE_max_Ha) bandE_max_Ha = WFK->eigen[k][n];
    }
  }
  
  /*convert E Ha min/max to eV and store*/
  ESTP->bandE_min = (bandE_min_Ha-fermi)*HATOEV;  
  ESTP->bandE_max = (bandE_max_Ha-fermi)*HATOEV;
  /*if want to cut off min and max at 10 eV, uncomment below*/
/*  if (ESTP->bandE_max > 10.0) {
    ESTP->bandE_max = 10.0;
  }
  if (ESTP->bandE_min < -10.0) {
    ESTP->bandE_min = -10.0;
  }
*/
  /*making min and max integers*/
  bandE_min_integer = floor(ESTP->bandE_min);
  bandE_max_integer = ceil(ESTP->bandE_max);
  
  /*find dE that corresponds to 0 eV or at the Ef*/
  zero = -(ESTP->bandE_min*EMESH)+0.5;
  ESTP->dE_zero = floor(zero);

  /*find the number of steps from Elow to Ehigh using energy mesh*/
  ESTP->nEstep = (EMESH*bandE_max_integer)-(EMESH*bandE_min_integer)+1;
  printf( "\tScanning energy range: %d eV to %d eV (%d steps)\n", bandE_min_integer, bandE_max_integer, ESTP->nEstep);
  printf( "\tdE step corresponding to 0 eV (E_F) = %d\n", ESTP->dE_zero);
 
} //END of find_energy_bounds function
  

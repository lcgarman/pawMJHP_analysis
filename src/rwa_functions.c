#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "rwa_functions.h"

void read_mjin_header(char filename[100], MottJonesConditions * MJC, FileCabinet* FCAB, FermiSphere* FS)
{
  /*This function reads in the COMPOUND.mjin file that contains information related to the mjhp analysis*/

  FILE* mjin_file; /*pointer to mjin file*/
  char strv[20]; /*arbitrary string*/
  char ABOfilename[100]; /*input filestring*/
  char MJOUTfilename[100]; /*output filestring*/
  int check; /*logic integer*/
  char lattice[10]; /*lattice and centering type*/
  double vec; /*number of valence electrons per prim uc*/
  int ndts; /*number of HKL datasets*/
  int n; /*index for dataset number*/
  int jzH, jzK, jzL; /*HKL indices of interest*/
  double scanE_min, scanE_mid, scanE_max; /*energy range of interest*/

  /*open the mjin file to be read*/
  printf("Reading %s....", filename);
  mjin_file=fopen(filename,"r+"); 
  if(mjin_file==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  /*first line is input file string*/
  fscanf(mjin_file, "%s\n", ABOfilename);
  strcpy(FCAB->ABOfilename, ABOfilename);
  /*second line is output file string*/
  fscanf(mjin_file, "%s\n", MJOUTfilename);
  strcpy(FCAB->MJOUTfilename, MJOUTfilename);
  /*third line is number of valence electrons in unit cell*/
  fscanf(mjin_file, "%s", strv);
  fscanf(mjin_file, "%lf", &vec);
  FS->vec = vec;
  /*fourth line is the lattice type and centering*/
  fscanf(mjin_file, "%s", lattice);
  strcpy(MJC->lattice, lattice);

  /*fifth line is number of HKL datasets*/
  fscanf(mjin_file, "%s", strv);
  fscanf(mjin_file, "%d", &ndts);
  MJC->ndts = ndts;
  
  /*check if this is the end of the file*/
  check=fscanf(mjin_file, "%s", strv);
  /*if so - finish reading and leave*/
  if (check==EOF) {
    fclose(mjin_file);
    return;
  }

  /*if not allocate and read in HKL indices*/
  MJC->jzH = AllocateMemory_oneD_int(MJC->jzH, ndts);
  MJC->jzK = AllocateMemory_oneD_int(MJC->jzK, ndts);
  MJC->jzL = AllocateMemory_oneD_int(MJC->jzL, ndts);
  for (n=0;n<ndts;n++) {
	fscanf(mjin_file, "%d", &jzH);
	fscanf(mjin_file, "%d", &jzK);
	fscanf(mjin_file, "%d", &jzL);
	MJC->jzH[n] = jzH;
	MJC->jzK[n] = jzK;
	MJC->jzL[n] = jzL;
  }
	
  /*again check if this is the end of the file*/
  check=fscanf(mjin_file, "%s", strv);
  /*if so - finish reading and leave*/
  if (check==EOF) {
    fclose(mjin_file);
    return;
  }

  /*if not allocate and read in HKL indices*/
  MJC->scanE_min = AllocateMemory_oneD_double(MJC->scanE_min, ndts);
  MJC->scanE_mid = AllocateMemory_oneD_double(MJC->scanE_mid, ndts);
  MJC->scanE_max = AllocateMemory_oneD_double(MJC->scanE_max, ndts);
  for (n=0;n<ndts;n++) {
	fscanf(mjin_file, "%lf", &scanE_min);
	fscanf(mjin_file, "%lf", &scanE_mid);
	fscanf(mjin_file, "%lf", &scanE_max);
	MJC->scanE_min[n] = scanE_min;
	MJC->scanE_mid[n] = scanE_mid;
	MJC->scanE_max[n] = scanE_max;
  }
	
  fclose(mjin_file);
  printf("Done.\n");
} //END of read_mjin_header 


void print_mjhpHKL_energy(char filename[200], VectorIndices *VECT, EnergyStep *ESTP, UnitCell *UC, EnergyContribution *ECON)
{
  /*This function prints out the mjhp HKL potential energy contributions*/
  FILE* fhkl; /*pointer to file*/
  int nEstep; /*number of energy steps*/
  int dE; /*energy step  index*/
  double Elow; /*lowest energy value*/
  int H, K, L; /*HKL indices*/
  double total_potential; /*total mjhp potential*/

  /*initialize variables*/
  nEstep = ESTP->nEstep;
  H = VECT->H;
  K = VECT->K;
  L = VECT->L;

  /*open file to print to*/
  fhkl = fopen(filename, "w");
  if(fhkl==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  printf("\nPrinting Energy for HKL to %s\n", filename);
  /*print header information*/
  fprintf(fhkl, "%s\n", filename);
  fprintf(fhkl, "astar: %lf\t%lf\t%lf\n", UC->ang_ax_star, UC->ang_ay_star, UC->ang_az_star);  
  fprintf(fhkl, "bstar: %lf\t%lf\t%lf\n", UC->ang_bx_star, UC->ang_by_star, UC->ang_bz_star);  
  fprintf(fhkl, "cstar: %lf\t%lf\t%lf\n", UC->ang_cx_star, UC->ang_cy_star, UC->ang_cz_star);  
  fprintf(fhkl, "noSteps %d \n", ESTP->nEstep);
  Elow = ESTP->bandE_min*EMESH;
  fprintf(fhkl, "lowestE %lf\n", floor(Elow));
  fprintf(fhkl, "Emesh %d\n", EMESH);
  fprintf(fhkl, "HKL %d\t%d\t%d\n", H, K, L);
  fprintf(fhkl, "\n");
  /*END of Header Info*/

  /*Print Energy by hkl*/
  total_potential = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    fprintf(fhkl, "%d\t%e\n", dE, ECON->total_potential[dE]);
    total_potential += ECON->total_potential[dE];
  }
  printf("\tTotal Potential Energy for %d %d %d = %e\n", VECT->H, VECT->K, VECT->L, total_potential);

  /*free energy contribution*/
  ECON->total_potential = FreeMemory_oneD_double(ECON->total_potential);

} //END of Print_Reflection function

void modify_line(char *line, Modification modifications[], int num_modifications)
{
  /*This function is doing the work of reading in one abinit in
 * file and printing a modified version to different abinit in file*/
  int max_line = 1024; /*maximum lenght of line to be read in*/
  int i; 
  char modified_line[max_line]; /*character containing a modified line*/

  /*search through number of modifications to be made*/
  for (i=0;i<num_modifications;i++) {
    /*search for the word that should be modified*/
    if (strstr(line, modifications[i].search) != NULL) {
      /*if found, make the modification - copy the line with the modification before*/
      snprintf(modified_line, sizeof(modified_line), "%s%s", modifications[i].replace, line);
      /*replace the original line with the modified version*/
      strcpy(line, modified_line);
      break;
    }
  } 

} /*END of modify_line*/
  
void modify_abinitin(char in_filename[100], char out_filename[100], TwoTheta * TTH)
{
  /*This function reads in the abinit in file and modifies it to perform
 * a non-selfconsistent ground state calculation at kpts of interest*/
  const int num_modifications = 8;  /*number of modifications to be made to file*/
  FILE* abin; /*abinit in file to read in*/
  FILE* mabin; /*abinit in file to be modified*/
  int max_line = 1024; /*largest line size that can be read in*/
  char line[max_line]; /*arbitrary character storing lines*/
  int nsym_bzk; /*number of symmetry-independent kpts*/
  int i; /*kpt index*/
  nsym_bzk = TTH->nsym_bzk;  

  /*initialize modifications structure*/
  Modification modifications[] = {
{"toldfe", "!"},
{"tolvrs", "!"},
{"occopt", "!"},
{"istwfk", "!"},
{"nshiftk", "!"},
{"shiftk", "!"},
{"ngkpt", "!"},
  };
  /*The above list will be modified in the new abinit in file 
 * by placing the ! before the word, so these lines are commented out*/

  /*open existing abinit in file to read*/
  printf("\nModifying %s\n", in_filename);
  abin=fopen(in_filename, "r");
  if (abin==NULL) {
    printf("%s not found. \n", in_filename);
    exit(0);
  }

  /*open a new file to write the modified in file too*/
  mabin=fopen(out_filename, "w");
  if (mabin==NULL) {
    printf("%s not found. \n", out_filename);
    exit(0);
  }

  /*read in the abinit in file and paste modified abinit in data to mabin*/
  while (fgets(line, sizeof(line), abin) != NULL) {
    modify_line(line, modifications, num_modifications);
    fputs(line, mabin);
  }

  /*add new kpt info lines to modified abinit in file*/
  fprintf(mabin,"\n!---- MJHP 2THETA VARS ----!\n");
  fprintf(mabin,"      istwfk %d*1\n", nsym_bzk);
  fprintf(mabin,"      tolwfr 1.0d-18\n");
  fprintf(mabin,"      nbdbuf -10\n");
  fprintf(mabin,"      iscf -2 \n");
  fprintf(mabin,"      kptopt 0\n");
  fprintf(mabin,"      nkpt %d\n", nsym_bzk);
  fprintf(mabin,"      kpt\n");
  /*print out kpts of interest*/
  for (i=0;i<nsym_bzk;i++) {
    fprintf(mabin,"      %lf %lf %lf\n", TTH->kx_sym[i], TTH->ky_sym[i], TTH->kz_sym[i]);
  }  
  /*free allocated TTH memory*/
  TTH->kx_sym = FreeMemory_oneD_double(TTH->kx_sym);
  TTH->ky_sym = FreeMemory_oneD_double(TTH->ky_sym);
  TTH->kz_sym = FreeMemory_oneD_double(TTH->kz_sym);

  /*close the original abinit in file and modified abinit in file*/
  fclose(abin); 
  fclose(mabin); 
      
} //END of  modify_abinitin function
    

void print_mjhp2theta_energy(char filename[200], TwoTheta *TTH, EnergyStep * ESTP, UnitCell* UC, EnergyContribution * ECON, FermiSphere *FS)
{
  /*This function prints the mjhp 2theta potential energy for each reflection*/
  FILE* f2th; /*pointer to mjout 2theta file*/
  int nEstep; /*number of energy bins*/
  int dE; /*energy bin index*/
  int nrflc; /*number of reflections*/
  int n; /*reflection index*/
  double Elow; /*minimum band energy*/
  double reflection_total; /*total summed energy from all reflections*/

  /*initialize variables*/
  nEstep = ESTP->nEstep;
  nrflc = TTH->nrflc;

  /*open the mjout file in write mode*/
  f2th = fopen(filename, "w");
  if(f2th==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  printf("\nPrinting energy for reflections to: %s\n", filename);
  /*print header information*/
  fprintf(f2th, "%s\n", filename);
  fprintf(f2th, "astar: %lf\t%lf\t%lf\n", UC->ang_ax_star, UC->ang_ay_star, UC->ang_az_star);  
  fprintf(f2th, "bstar: %lf\t%lf\t%lf\n", UC->ang_bx_star, UC->ang_by_star, UC->ang_bz_star);  
  fprintf(f2th, "cstar: %lf\t%lf\t%lf\n", UC->ang_cx_star, UC->ang_cy_star, UC->ang_cz_star);  
  fprintf(f2th, "noSteps %d \n", ESTP->nEstep);
  Elow = ESTP->bandE_min*EMESH;
  fprintf(f2th, "lowestE %lf\n", floor(Elow));
  fprintf(f2th, "Emesh %d\n", EMESH);
  fprintf(f2th, "FS_angle %lf\n", FS->two_theta);
  fprintf(f2th, "nrflc %d\n", nrflc);
  fprintf(f2th, "\n");
  /*END of Header Info*/

  /*Print energy by reflection*/
  for (n=0;n<nrflc;n++) {
    fprintf(f2th, "\n");
    fprintf(f2th, "rflc#%d\t %lf\t %d %d %d\t%d\n", n, TTH->two_theta[n], TTH->rflc_H[n], TTH->rflc_K[n], TTH->rflc_L[n], TTH->rflc_mult[n]);
    reflection_total = 0.0;
    for (dE=0;dE<nEstep;dE++) {
      fprintf(f2th, "%d\t%e\n", dE, ECON->rflc_total[dE][n]);
      reflection_total += ECON->rflc_total[dE][n];
    }
    printf("\tTotal Potential Energy for %d %d %d = %e\n", TTH->rflc_H[n], TTH->rflc_K[n], TTH->rflc_L[n], reflection_total);
  }

  /*close file*/
  fclose(f2th); 
  /*free energy contribution*/
  ECON->rflc_total = FreeMemory_twoD_double(ECON->rflc_total, nEstep);

} //END of print_mjhp2theta_energy function


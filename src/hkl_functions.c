#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "hkl_functions.h"

void transform_HKL(MottJonesConditions * MJC) 
{
  /*This function transforms the HKL indices from conventional to primitive.*/
  char lattice[10]; /*bravais lattice type*/
  int len_lattice; /*lenght of lattice varibale above*/
  char center; /*centering type*/
  /*characters for each of the centerings*/
  char primitive = 'P';
  char bodycenter = 'I';
  char facecenter = 'F';
  char acenter = 'A';
  char bcenter = 'B';
  char ccenter = 'C';
  int H_symm, K_symm, L_symm; /*HKL indices after the conversion to primitive*/
  int H, K, L; /*varibale that stores primitive indices*/
  
  /*storing conventional HKL indices*/
  H = MJC->nH; 
  K = MJC->nK;
  L = MJC->nL;

  /*separating lattice variable. ie - "cP" -> "c" and "P"*/
  strcpy(lattice, MJC->lattice);
  len_lattice = strlen(lattice);
  /*check lattice is right length*/
  if (len_lattice > 2) printf("ERROR: length of lattice (%d) is too long.\n", len_lattice);
  center = lattice[1];

  printf("Transforming HKL indices to primitive centering.\n");
  /*Transform to primitive cell*/
  if (center==primitive) {
    printf("\tAlready in primitive centering.\n");
  }
  /*transormation for body centering*/
  else if (center==bodycenter) { 
    printf("\nTransforming HKL: I->P\n");
	H_symm = (-0.5*H)+(0.5*K)+(0.5*L);
	K_symm = (0.5*H)+(-0.5*K)+(0.5*L);
	L_symm = (0.5*H)+(0.5*K)+(-0.5*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transormation for face centering*/
  else if (center==facecenter) { 
    printf("\nTransforming HKL: F->P\n");
	H_symm = (0.0*H)+(0.5*K)+(0.5*L);
	K_symm = (0.5*H)+(0.0*K)+(0.5*L);
	L_symm = (0.5*H)+(0.5*K)+(0.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transormation for a-centering*/
  else if (center==acenter) { 
    printf("\nTransforming HKL: A->P\n");
	H_symm = (1.0*H)+(0.0*K)+(0.0*L);
	K_symm = (0.0*H)+(0.5*K)+(-0.5*L);
	L_symm = (0.0*H)+(0.5*K)+(0.5*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transormation for b-centering*/
  else if (center==bcenter) { 
    printf("\nTransforming HKL: B->P\n");
	H_symm = (0.5*H)+(0.0*K)+(-0.5*L);
	K_symm = (0.0*H)+(1.0*K)+(0.0*L);
	L_symm = (0.5*H)+(0.0*K)+(0.5*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transormation for c-centering*/
  else if (center==ccenter) { 
    printf("\nTransforming HKL: C->P\n");
	H_symm = (0.5*H)+(-0.5*K)+(0.0*L);
	K_symm = (0.5*H)+(0.5*K)+(0.0*L);
	L_symm = (0.0*H)+(0.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else {
    printf("ERROR: No Code written for %s\n", MJC->lattice);
    printf("\t case sensitive (xY) \n");
    exit(0);
  }

  /*store primitive indices as the HKL indices*/
  MJC->nH = H;
  MJC->nK = K;
  MJC->nL = L;

} //END of Transform_hklreflection

void find_symmetric_hkl(VectorIndices* VECT, Symmetry* SYM, NumberGrid* GRD) 
{
  /*This function finds all symmetry equivalent HKL indices.*/
  int j; 
  int sym; /*index for symmetry*/
  int nsym; /*number of symmetry elements*/
  int H, K, L; /*HKL indices*/
  int nHKL; /*number of symmetry equivalent indices*/
  int nomatch; /*integer to track number of if matches are found*/
  int H_symm, K_symm, L_symm; /*symmetrized HKL indices*/
  int* H_arr; /*array of symmetric HKL indices*/
  int* K_arr;
  int* L_arr;
  int ngfftx, ngffty, ngfftz; /*number grid*/
  
  /*initializing variables*/
  H = VECT->H;
  K = VECT->K;
  L = VECT->L;
  nsym = SYM->nsym;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  H_arr = NULL;
  K_arr = NULL;
  L_arr = NULL;
  
  /*Allocate Memory for HKL arrays*/
  H_arr = AllocateMemory_oneD_int(H_arr, nsym);
  K_arr = AllocateMemory_oneD_int(K_arr, nsym);
  L_arr = AllocateMemory_oneD_int(L_arr, nsym);

  printf( "\nSymmetrically equivalent HKL:\n");
  /*assign first element of HKL array to input HKL*/
  H_arr[0] = H;
  K_arr[0] = K;
  L_arr[0] = L;
  nHKL = 1;
  /*Find symmetry equivalent reflections based off of nonsymmorphic symm elements*/
  for(sym=0;sym<nsym;sym++) {
	H_symm = SYM->symrel[0][0][sym]*H + SYM->symrel[1][0][sym]*K + SYM->symrel[2][0][sym]*L;
	K_symm = SYM->symrel[0][1][sym]*H + SYM->symrel[1][1][sym]*K + SYM->symrel[2][1][sym]*L;
	L_symm = SYM->symrel[0][2][sym]*H + SYM->symrel[1][2][sym]*K + SYM->symrel[2][2][sym]*L;

	/*search elements in array if you find a match go to next symm element; if not store*/
	nomatch = 0;
    for (j=0;j<nHKL;j++) {
	  if ((H_symm == H_arr[j])&&(K_symm == K_arr[j])&&(L_symm == L_arr[j])) {
	    break;
	  }
	  else nomatch++;
    }
	if (nomatch == nHKL) {
	  H_arr[nHKL] = H_symm;
	  K_arr[nHKL] = K_symm;
	  L_arr[nHKL] = L_symm;
	  nHKL++;
	}
  } /*end of sym->nsym loop*/
  
  printf( "\tMultiplicity of HKL Reflection = %d\n", nHKL);
  /*allocate memory for RFLC arrays*/
  VECT->H_arr = AllocateMemory_oneD_int(VECT->H_arr, nHKL);
  VECT->K_arr = AllocateMemory_oneD_int(VECT->K_arr, nHKL);
  VECT->L_arr = AllocateMemory_oneD_int(VECT->L_arr, nHKL);
  VECT->H_posarr = AllocateMemory_oneD_int(VECT->H_posarr, nHKL);
  VECT->K_posarr = AllocateMemory_oneD_int(VECT->K_posarr, nHKL);
  VECT->L_posarr = AllocateMemory_oneD_int(VECT->L_posarr, nHKL);
  /*end of allocation*/

  /*store variables in VectorIndices structure*/
  for (j=0;j<nHKL;j++) {
    VECT->H_arr[j] = H_arr[j];	
    if (H_arr[j] < 0) {
      VECT->H_posarr[j] = H_arr[j] + ngfftx; 
    }
    else {
      VECT->H_posarr[j] = H_arr[j]; 
    }
    VECT->K_arr[j] = K_arr[j]; 
    if (K_arr[j] < 0) {
      VECT->K_posarr[j] = K_arr[j] + ngffty; 
    }
    else {
      VECT->K_posarr[j] = K_arr[j]; 
    }
    VECT->L_arr[j] = L_arr[j];
    if (L_arr[j] < 0) {
      VECT->L_posarr[j] = L_arr[j] + ngfftz; 
    }
    else {
      VECT->L_posarr[j] = L_arr[j]; 
    }
	printf( "\t#%d:   %d %d %d \t (%d %d %d)\n", j+1, VECT->H_arr[j], VECT->K_arr[j], VECT->L_arr[j], VECT->H_posarr[j], VECT->K_posarr[j], VECT->L_posarr[j]);
  }
  VECT->nHKL = nHKL;

  /*free temporary allocated variables*/
  H_arr = FreeMemory_oneD_int(H_arr);
  K_arr = FreeMemory_oneD_int(K_arr);
  L_arr = FreeMemory_oneD_int(L_arr);

} //END of Find_HKLsymmetry

Conventional HKL_convert_toP(MottJonesConditions * MJC, int H, int K, int L)
{
  /*This function returns the conventional indices 
 * from a structure contaning primitve indices.*/
  char lattice[10]; /*bravais lattice type*/
  int len_lattice; /*length of lattice name*/
  char center; /*centering of lattice*/
  char primitive = 'P';  /*different types of centering*/
  char bodycenter = 'I';
  char facecenter = 'F';
  char acenter = 'A';
  char bcenter = 'B';
  char ccenter = 'C';
  int H_symm, K_symm, L_symm; /*symmetric HKL indices based off centering conversion*/
  
  /*initializing structure*/
  Conventional con;

  /*separating lattice variable. ie - "cP" -> "c" and "P"*/
  strcpy(lattice, MJC->lattice);
  len_lattice = strlen(lattice);

  /*check lattice is right length*/
  if (len_lattice > 2) printf("ERROR: length of lattice (%d) is too long.\n", len_lattice);
  center = lattice[1];

  /*Transform to primitive cell*/
  if (center==primitive) {
    printf("Already in Conventional Centering (P)\n");
  }
  /*transform for body-centering*/
  else if (center==bodycenter) { 
    printf("Transforming HKL: P->I\n");
	H_symm = (0.0*H)+(1.0*K)+(1.0*L);
	K_symm = (1.0*H)+(0.0*K)+(1.0*L);
	L_symm = (1.0*H)+(1.0*K)+(0.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transform for face-centering*/
  else if (center==facecenter) { 
    printf("Transforming HKL: P->F\n");
	H_symm = (-1.0*H)+(1.0*K)+(1.0*L);
	K_symm = (1.0*H)+(-1.0*K)+(1.0*L);
	L_symm = (1.0*H)+(1.0*K)+(-1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transform for a-centering*/
  else if (center==acenter) { 
    printf("Transforming HKL: A->P\n");
	H_symm = (1.0*H)+(0.0*K)+(0.0*L);
	K_symm = (0.0*H)+(1.0*K)+(1.0*L);
	L_symm = (0.0*H)+(-1.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transform for b-centering*/
  else if (center==bcenter) { 
    printf("Transforming HKL: P->B\n");
	H_symm = (1.0*H)+(0.0*K)+(1.0*L);
	K_symm = (0.0*H)+(1.0*K)+(0.0*L);
	L_symm = (-1.0*H)+(0.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  /*transform for c-centering*/
  else if (center==ccenter) { 
    printf("Transforming HKL: P->C\n");
	H_symm = (1.0*H)+(1.0*K)+(0.0*L);
	K_symm = (-1.0*H)+(1.0*K)+(0.0*L);
	L_symm = (0.0*H)+(0.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else {
    printf("ERROR: No Code written for %s\n", MJC->lattice);
    printf("\t case sensitive (xY) \n");
    exit(0);
  }

  /*store conventional indices*/
  con.H = H;
  con.K = K;
  con.L = L;
  
  /*return the conventional indices*/
  return(con);

} //END of primitive_to_conventional 


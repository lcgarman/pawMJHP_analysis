#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structures.h"
#include "allocate_memory.h"

void Initialize_NumberGrid(NumberGrid * NG) 
{
  NG->h_grid = NULL;
  NG->k_grid = NULL;
  NG->l_grid = NULL;
}

void Initialize_Symmetry(Symmetry * SYM) 
{
  SYM->symrel = NULL;
  SYM->mult = NULL;
  SYM->tnons = NULL;
}
 
void Initialize_Wavefunction(Wavefunction * WFK)
{
  WFK->npw = NULL;
  WFK->kpt = NULL;
  WFK->wtk = NULL;
  WFK->kg = NULL;
  WFK->eigen = NULL;
  WFK->occ = NULL;
  WFK->cg = NULL;
}

void Initialize_BinaryGrid(BinaryGrid * BIN)
{
  BIN->real_grid = NULL;
  BIN->rec_grid = NULL;
  BIN->cc_rec_grid = NULL;
}

void Initialize_AtomicVariables(AtomicVariables * ATM) 
{
  ATM->xred = NULL;
  ATM->typat = NULL;
  ATM->atomicno = NULL;
  ATM->Xcart = NULL;
  ATM->Ycart = NULL;
  ATM->Zcart = NULL;
}

void Initialize_VectorIndices(VectorIndices * VECT) 
{
  VECT->H_arr = NULL;
  VECT->K_arr = NULL;
  VECT->L_arr = NULL;
  VECT->H_posarr = NULL;
  VECT->K_posarr = NULL;
  VECT->L_posarr = NULL;
}

void Initialize_MottJonesConditions(MottJonesConditions * MJC)
{
  MJC->jzH = NULL;
  MJC->jzK = NULL;
  MJC->jzL = NULL;
  MJC->JZkpt = NULL;
  MJC->JZk_sym = NULL;
  MJC->JZmult = NULL;
  MJC->JZwtk = NULL;
  MJC->BZk = NULL;
  MJC->Hpw = NULL;
  MJC->Kpw = NULL;
  MJC->Lpw = NULL;
  MJC->scanE_min = NULL;
  MJC->scanE_mid = NULL;
  MJC->scanE_max = NULL;
}

void Initialize_TwoTheta(TwoTheta * TTH)
{
  TTH->two_theta = NULL;
  TTH->d_hkl = NULL;
  TTH->F_hkl = NULL;
  TTH->rflc_H = NULL;
  TTH->rflc_K = NULL;
  TTH->rflc_L = NULL;
  TTH->hpw = NULL;
  TTH->kpw = NULL;
  TTH->lpw = NULL;
  TTH->rflc_mult = NULL;
  TTH->kx = NULL;
  TTH->ky = NULL;
  TTH->kz = NULL;
  TTH->kx_sym = NULL;
  TTH->ky_sym = NULL;
  TTH->kz_sym = NULL;
  TTH->BZkpt = NULL;
  TTH->BZkpt_sym = NULL;
  TTH->hpw_sym = NULL;
  TTH->kpw_sym = NULL;
  TTH->lpw_sym = NULL;
  TTH->rflc_H_sym = NULL;
  TTH->rflc_K_sym = NULL;
  TTH->rflc_L_sym = NULL;
}

void Initialize_EnergyContribution(EnergyContribution * ECON)
{
  ECON->local = NULL;
  ECON->nonlocal = NULL;
  ECON->kinetic = NULL;
  ECON->KE_correction = NULL;
  ECON->total_potential = NULL;
  ECON->rflc_local = NULL;
  ECON->rflc_nonlocal = NULL;
  ECON->rflc_KE_correction = NULL;
  ECON->rflc_total = NULL;
}

void Initialize_PawAtomicData(PawAtomicData * PAW) 
{
  PAW->nrhoij = NULL;
  PAW->irhoij = NULL;
  PAW->rhoij = NULL;
  PAW->Dij = NULL;
  PAW->Qhat = NULL;
  PAW->lnmax = NULL;
  PAW->lmn_size = NULL;
  PAW->mmax = NULL;
  PAW->psp_path = NULL;
  PAW->l_orbital =  NULL;
  PAW->paw_rad = NULL;
  PAW->shape_fxn_rc = NULL;
  PAW->paw_lnrad = NULL;
  PAW->values = NULL;
  PAW->derivatives = NULL;
  PAW->ae_core_density = NULL;
  PAW->pseudo_core_density = NULL;
  PAW->pseudo_valence_density = NULL;
  PAW->zero_potential = NULL;
  PAW->local_ionic_potential = NULL;
  PAW->ae_partial_wave = NULL;
  PAW->pseudo_partial_wave = NULL;
  PAW->projector_function = NULL;
  PAW->kinetic_energy_differences = NULL;
  PAW->KE_diff_matrix = NULL;
  PAW->exact_exchange_matrix = NULL;
  PAW->Dij_matrix = NULL;
  PAW->c_squig = NULL;
}

int* AllocateMemory_oneD_int(int *array, int dim1)
{
  /*allocates memory for a 1D arry of integers*/
  int i;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }
    
  /*allocate 1D array to pointer*/
  array = malloc(dim1 * sizeof(int));
  /*check allocation was successful*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) array[i] = 0;
  /*return the array at end of function*/
  return array;
} 

int** AllocateMemory_twoD_int(int **array, int dim1, int dim2)
{
  /*allocates memory for a 2D arry of integers*/
  int i, j;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int*));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed (1D/2)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(int));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed(2D/2)\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      array[i][j] = 0;
    }
  }

  /*return the array at end of function*/
  return array;
} 

int*** AllocateMemory_threeD_int(int ***array, int dim1, int dim2, int dim3)
{
  /*allocates memory for a 3D arry of integers*/
  int i, j, k;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(int*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(int));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        array[i][j][k] = 0;
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 

int**** AllocateMemory_fourD_int(int ****array, int dim1, int dim2, int dim3, int dim4)
{
  /*allocates memory for a 4D arry of integers*/
  int i, j, k, l;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int***));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(int**));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(int*));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }

      /*allocate fourth dimension*/
      for (k=0;k<dim3;k++) {
        array[i][j][k] = malloc(dim4 * sizeof(int));
        /*check allocation*/
		if (array[i][j][k] == NULL) {
		  printf("ERROR: Memory Allocation Failed\n");
		  for (l=0;l<k;l++) {
			free(array[i][j][l]);
		  }
		  free(array[i][j]);
		  free(array[i]);
		  free(array);
		  exit(0);
		}
	  }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        for (l=0;l<dim4;l++) {
          array[i][j][k][l] = 0;
        }
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 

double* AllocateMemory_oneD_double(double *array, int dim1)
{
  /*allocates memory for a 1D arry of doubles*/
  int i;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
    printf("NOT NULL\n");
  }

  /*allocate 1D array to pointer*/
  array = malloc(dim1 * sizeof(double));
  /*check allocation was successful*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
	array[i] = 0;
  }

  /*return the array at end of function*/
  return array;
} 

double** AllocateMemory_twoD_double(double **array, int dim1, int dim2)
{
  /*allocates memory for a 2D arry of doubles*/
  int i, j;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of array*/
  array = malloc(dim1 * sizeof(double*));
  /*check allocation for every dimension*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed (1D)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(double));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
	  array[i][j] = 0.0;
    }
  }

  return array;
}

double*** AllocateMemory_threeD_double(double ***array, int dim1, int dim2, int dim3)
{
  /*allocates memory for a 3D arry of doubles*/
  int i, j, k;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }


  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(double*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(double));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
		array[i][j][k] = 0.0;
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 

double**** AllocateMemory_fourD_double(double ****array, int dim1, int dim2, int dim3, int dim4)
{
  /*allocates memory for a 4D arry of doubles*/
  int i, j, k, l;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double***));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(double**));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(double*));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
      
      /*allocate fourth dimension*/
      for (k=0;k<dim3;k++) {
        array[i][j][k] = malloc(dim4 * sizeof(double));
        /*check allocation*/
        if (array[i][j][k] == NULL) {
          printf("ERROR: Memory Allocation Failed\n");
          for (l=0;l<k;l++) {
            free(array[i][j][l]);
          }
          free(array[i][j]);
          free(array[i]);
          free(array);
          exit(0);
        }
	  }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        for (l=0;l<dim4;l++) {
		  array[i][j][k][l] = 0.0;
        }
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 

gsl_complex*** AllocateMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2, int dim3)
{
  /*allocates memory for a 3D arry of complex doubles*/
  int i, j, k;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(gsl_complex**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(gsl_complex*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(gsl_complex));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        GSL_REAL(array[i][j][k]) = 0.0;
        GSL_IMAG(array[i][j][k]) = 0.0;
      }
    }
  }
        
  /*return the array at end of function*/
  return array;
} 

gsl_complex**** AllocateMemory_fourD_complex(gsl_complex**** array, int dim1, int dim2, int dim3, int dim4)
{
  /*allocates memory for a 4D arry of doubles*/
  int i, j, k, l;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  printf("Allocating for 4D complex variable...");
  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(gsl_complex***));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(gsl_complex**));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(gsl_complex*));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }

      /*allocate fourth dimension*/
      for (k=0;k<dim3;k++) {
        array[i][j][k] = malloc(dim4 * sizeof(gsl_complex));
        /*check allocation*/
        if (array[i][j][k] == NULL) {
          printf("ERROR: Memory Allocation Failed\n");
          for (l=0;l<k;l++) {
            free(array[i][j][l]);
          }
          free(array[i][j]);
          free(array[i]);
          free(array);
          exit(0);
        }
	  }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        for (l=0;l<dim4;l++) {
		  GSL_REAL(array[i][j][k][l]) = 0.0;
		  GSL_IMAG(array[i][j][k][l]) = 0.0;
        }
      }
    }
  }

  /*return the array at end of function*/
  printf(" successful.\n");
  return array;
}

gsl_complex***** AllocateMemory_fiveD_complex(gsl_complex***** array, int dim1, int dim2, int dim3, int dim4, int dim5)
{
  /*allocates memory for a 5D arry of complex doubles*/
  int i, j, k, l, m;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(gsl_complex****));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(gsl_complex***));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(gsl_complex**));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }

      /*allocate fourth dimension*/
      for (k=0;k<dim3;k++) {
        array[i][j][k] = malloc(dim4 * sizeof(gsl_complex*));
        /*check allocation*/
        if (array[i][j][k] == NULL) {
          printf("ERROR: Memory Allocation Failed\n");
          for (l=0;l<k;l++) {
            free(array[i][j][l]);
          }
          free(array[i][j]);
          free(array[i]);
          free(array);
          exit(0);
        }
  
        /*allocate fifth dimension*/
        for (l=0;l<dim4;l++) {
          array[i][j][k][l] = malloc(dim5 * sizeof(gsl_complex));
          /*check allocation*/
          if (array[i][j][k][l] == NULL) {
            for (m=0;m<l;m++) {
              free(array[i][j][k][m]);
            }
            free(array[i][j][k]);
            free(array[i][j]);
            free(array[i]);
            free(array);
            exit(0);
          }
        }
	  }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        for (l=0;l<dim4;l++) {
          for (m=0;m<dim5;m++) {
     	    GSL_REAL(array[i][j][k][l][m]) = 0.0;
     	    GSL_IMAG(array[i][j][k][l][m]) = 0.0;
          }
        }
      }
    }
  }

  /*return the array at end of function*/
  return array;
}


int* FreeMemory_oneD_int(int *array)
{
  /*frees memory from a 1D arry of integers*/

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  free(array);
  array = NULL;
  return array;
}

int** FreeMemory_twoD_int(int **array, int dim1)
{
  /*frees memory from a 2D arry of integers*/
  int i;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

int*** FreeMemory_threeD_int(int ***array, int dim1, int dim2)
{
  /*frees memory from a 3D arry of integers*/
  int i, j;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

int**** FreeMemory_fourD_int(int ****array, int dim1, int dim2, int dim3)
{
  /*frees memory from a 4D arry of integers*/
  int i, j, k;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        free(array[i][j][k]);
      }
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

double* FreeMemory_oneD_double(double *array)
{
  /*frees memory from a 1D arry of doubles*/

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  free(array);
  array = NULL;
  return array;
}

double** FreeMemory_twoD_double(double **array, int dim1)
{
  /*frees memory from a 2D arry of doubles*/
  int i;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

double*** FreeMemory_threeD_double(double ***array, int dim1, int dim2)
{
  /*frees memory from a 3D arry of doubles*/
  int i, j;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;

  return array;
}

double**** FreeMemory_fourD_double(double ****array, int dim1, int dim2, int dim3)
{
  /*frees memory from a 4D arry of doubles*/
  int i, j, k;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        free(array[i][j][k]);
      }
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}


gsl_complex*** FreeMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2)
{
  /*frees memory from a 3D arry of complex doubles*/
  int i, j;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

gsl_complex**** FreeMemory_fourD_complex(gsl_complex ****array, int dim1, int dim2, int dim3)
{
  /*frees memory from a 4D arry of complex doubles*/
  int i, j, k;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
	    free(array[i][j][k]);
	  }
	  free(array[i][j]);
	}
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

gsl_complex***** FreeMemory_fiveD_complex(gsl_complex *****array, int dim1, int dim2, int dim3, int dim4)
{
  /*frees memory from a 5D arry of complex doubles*/
  int i, j, k, l;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        for (l=0;l<dim4;l++) {
          free(array[i][j][k][l]);
		} 
	    free(array[i][j][k]);
	  }
	  free(array[i][j]);
	}
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

char** AllocateMemory_twoD_char(char **array, int dim1, int dim2)
{
  /*allocates memory from a 2D arry of characters*/
  int i, j;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(char*));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed (1D/2)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(char));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed(2D/2)\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*NO ZERO OUT of MEMORY for CHARACTERS IMPLEMENTED*/  
  /*return the array at end of function*/
  return array;
}

char** FreeMemory_twoD_char(char **array, int dim1)
{
  /*frees memory from a 2D arry of characters*/
  int i;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    free(array[i]);
  }
  free(array);

  return array;
}


void AllocateMemory_Wavefunctions(Wavefunction* WFK) 
{
  /*allocates memory for wavefunction coefficients and coordinates file*/
  int nkpt;
  int nband;
  int npw;
  int i, j, k, l;

  nkpt = WFK->nkpt;
  nband = WFK->nband;

  if (WFK->cg != NULL) {
    printf("ERROR: wavefunction coefficients not empty before allocation\n");
  }
  if (WFK->kg != NULL) {
    printf("ERROR: wavefunction coordinates not empty before allocation\n");
  }

  /*Allocating Memory for Wavefunction Coefficients cg*/
  WFK->cg = malloc(nkpt * sizeof(double***));
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    WFK->cg[i] = malloc(nband * sizeof(double**));
    for (j=0;j<nband;j++) {
      WFK->cg[i][j] = malloc(npw * sizeof(double*));
	  for (k=0;k<npw;k++) {
        WFK->cg[i][j][k] = malloc(2 * sizeof(double));
	  }
	}
  }
  /*Checking if Memory was allocated successfully*/
  if (WFK->cg == NULL) {
	printf("ERROR: Memory Allocation Failed for Wavefunction Coefficients (cg) \n");
    exit(0);
  } 
  
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    for (j=0;j<nband;j++) {
      for (k=0;k<npw;k++) {
        for (l=0;l<2;l++) {
          WFK->cg[i][j][k][l] = 0.0;
		}
	  }
	}
  }

  /*Allocating Memory for Reduced plane wave Coordinates (kg)*/
  WFK->kg = malloc(nkpt * sizeof(int**));
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    WFK->kg[i] = malloc(npw * sizeof(int*));
    for (k=0;k<npw;k++) {
      WFK->kg[i][k] = malloc(3 * sizeof(int));
	}
  }
  if (WFK->kg == NULL) {
	printf("\nERROR: Memory Allocation Failed for PW Coordinates (kg) \n");
    exit(0);
  }
  
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    for (k=0;k<npw;k++) {
      for (l=0;l<3;l++) {
        WFK->kg[i][k][l] = 0.0;
      }
    }
  }
  
}
void FreeMemory_Wavefunctions(Wavefunction* WFK) 
{
  /*frees memory related to wavefunction*/
  int nkpt;
  int nband;
  int npw;
  int i, j, k;

  nkpt = WFK->nkpt;
  nband = WFK->nband;

  if ((WFK->kg == NULL)) {
    return;
  }
  if ((WFK->cg == NULL)) {
    return;
  }

  /*free planewave reduced coordinates*/
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    for (k=0;k<npw;k++) {
      free(WFK->kg[i][k]);
    }
    free(WFK->kg[i]);
  }
  free(WFK->kg);
  WFK->kg = NULL;

  /*free wavefunction coefficients*/
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    for (j=0;j<nband;j++) {
      for (k=0;k<npw;k++) {
        free(WFK->cg[i][j][k]);
      }
      free(WFK->cg[i][j]);
    }
    free(WFK->cg[i]);
  }
  free(WFK->cg);
  WFK->cg = NULL;

}

void AllocateMemory_PAWvariable(PawAtomicData* PAW, AtomicVariables * ATM)
{
  /*function that controls the allocation of the PAW related variables*/

  //OPTION==0 -> mmax as second variable
  //OPTION==1 -> lmn_size as second variable

  /*allocate memory for ALL paw variables*/
  PAW->paw_rad = AllocateMemory_oneD_double(PAW->paw_rad, ATM->ntypat);
  PAW->shape_fxn_rc = AllocateMemory_oneD_double(PAW->shape_fxn_rc, ATM->ntypat);
  PAW->paw_lnrad = AllocateMemory_PAWvariable_twoD(PAW->paw_lnrad, PAW, ATM->ntypat, 1);
  PAW->values =  AllocateMemory_PAWvariable_twoD(PAW->values, PAW, ATM->ntypat, 0);
  PAW->derivatives = AllocateMemory_PAWvariable_twoD(PAW->derivatives, PAW, ATM->ntypat, 0);
  PAW->ae_core_density = AllocateMemory_PAWvariable_twoD(PAW->ae_core_density, PAW, ATM->ntypat, 0);
  PAW->pseudo_core_density = AllocateMemory_PAWvariable_twoD(PAW->pseudo_core_density, PAW, ATM->ntypat, 0);
  PAW->pseudo_valence_density = AllocateMemory_PAWvariable_twoD(PAW->pseudo_valence_density, PAW, ATM->ntypat, 0);
  PAW->zero_potential = AllocateMemory_PAWvariable_twoD(PAW->zero_potential, PAW, ATM->ntypat, 0);
  PAW->local_ionic_potential = AllocateMemory_PAWvariable_twoD(PAW->local_ionic_potential, PAW, ATM->ntypat, 0);
  PAW->kinetic_energy_differences = AllocateMemory_PAWvariable_twoD(PAW->kinetic_energy_differences, PAW, ATM->ntypat, 2);
  PAW->exact_exchange_matrix = AllocateMemory_PAWvariable_twoD(PAW->exact_exchange_matrix, PAW, ATM->ntypat, 2);
  PAW->ae_partial_wave = AllocateMemory_PAWvariable_threeD(PAW->ae_partial_wave, PAW, ATM->ntypat);
  PAW->pseudo_partial_wave = AllocateMemory_PAWvariable_threeD(PAW->pseudo_partial_wave, PAW, ATM->ntypat);
  PAW->projector_function = AllocateMemory_PAWvariable_threeD(PAW->projector_function, PAW, ATM->ntypat);

  /*requires different memory allocation type*/
  PAW->KE_diff_matrix = AllocateMemory_PAWthreeD_double_type2(PAW->KE_diff_matrix, PAW, ATM); 

  printf("Memory Allocated successfully for PAW variables.\n");
} 

void FreeMemory_PAWvariable(PawAtomicData* PAW, AtomicVariables * ATM)
{
  /*Frees memory for PAWvariables*/

  /*Free memory for 3D Dij and Rhoij matrices*/
  PAW->Dij_matrix = FreeMemory_PAWthreeD_double(PAW->Dij_matrix, PAW, ATM);
  PAW->KE_diff_matrix = FreeMemory_PAWthreeD_double_type2(PAW->KE_diff_matrix, PAW, ATM); 
  
  /* Free memory for twoD arrays*/
  /*do not need specific function as second dimension is not needed*/
  PAW->values = FreeMemory_twoD_double(PAW->values, ATM->ntypat); 
  PAW->derivatives = FreeMemory_twoD_double(PAW->derivatives, ATM->ntypat); 
  PAW->ae_core_density = FreeMemory_twoD_double(PAW->ae_core_density, ATM->ntypat); 
  PAW->pseudo_core_density = FreeMemory_twoD_double(PAW->pseudo_core_density, ATM->ntypat); 
  PAW->pseudo_valence_density = FreeMemory_twoD_double(PAW->pseudo_valence_density, ATM->ntypat); 
  PAW->zero_potential = FreeMemory_twoD_double(PAW->zero_potential, ATM->ntypat); 
  PAW->local_ionic_potential = FreeMemory_twoD_double(PAW->local_ionic_potential, ATM->ntypat); 
  PAW->kinetic_energy_differences = FreeMemory_twoD_double(PAW->kinetic_energy_differences, ATM->ntypat); 
  PAW->exact_exchange_matrix = FreeMemory_twoD_double(PAW->exact_exchange_matrix, ATM->ntypat); 

  /*Free mememory for threeD arrays*/
  PAW->ae_partial_wave = FreeMemory_PAWvariable_threeD(PAW->ae_partial_wave, PAW, ATM->ntypat);
  PAW->pseudo_partial_wave = FreeMemory_PAWvariable_threeD(PAW->ae_partial_wave, PAW, ATM->ntypat);
  PAW->projector_function = FreeMemory_PAWvariable_threeD(PAW->ae_partial_wave, PAW, ATM->ntypat);

  /*free memory for oneD arrays*/
  PAW->lnmax = FreeMemory_oneD_int(PAW->lnmax);
  PAW->lmn_size = FreeMemory_oneD_int(PAW->lmn_size);
  PAW->l_orbital = FreeMemory_twoD_int(PAW->l_orbital, ATM->ntypat);
  PAW->mmax = FreeMemory_oneD_int(PAW->mmax);
  PAW->psp_path = FreeMemory_twoD_char(PAW->psp_path, ATM->ntypat);
  PAW->paw_lnrad = FreeMemory_twoD_double(PAW->paw_lnrad, ATM->ntypat); 
  PAW->paw_rad = FreeMemory_oneD_double(PAW->paw_rad);
  PAW->shape_fxn_rc = FreeMemory_oneD_double(PAW->shape_fxn_rc);
  
  /*FreeMemory for paw matrices rhoij and Dij*/
  PAW->nrhoij = FreeMemory_oneD_int(PAW->nrhoij);
  PAW->irhoij = FreeMemory_twoD_int(PAW->irhoij, ATM->natom);
  PAW->rhoij = FreeMemory_twoD_double(PAW->rhoij, ATM->natom);
  PAW->Dij = FreeMemory_twoD_double(PAW->Dij, ATM->natom);

}

double** AllocateMemory_PAWvariable_twoD(double** array, PawAtomicData* PAW, int dim1, int option)
{
  /*allocates 2D array of double specific to PAW atomic dataset;
 * 1st dimension is given as input (usually ntypat)
 * 2nd dimension option=0->mmax and option=1->lmn_size and option=2->lnmax*/

  int dim2;
  int i;
  int j;
  
  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of array*/
  array = malloc(dim1 * sizeof(double*));
  /*check allocation for every dimension*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed (1D)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    if (option==0) dim2 = PAW->mmax[i];
    else if (option==2) dim2 = PAW->lmn_size[i]*2;
    else if (option==1) dim2 = PAW->lnmax[i];
    array[i] = malloc(dim2 * sizeof(double));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    if (option==0) dim2 = PAW->mmax[i];
    else if (option==1) dim2 = PAW->lnmax[i];
    else if (option==2) dim2 = PAW->lmn_size[i];
    for (j=0;j<dim2;j++) {
	  array[i][j] = 0.0;
    }
  }

  return array;
}

double*** AllocateMemory_PAWvariable_threeD(double*** array, PawAtomicData* PAW, int dim1)
{
  /*allocates 3D array of double specific to PAW atomic dataset;
 * 1st dimension is given as input (usually ntypat)
 * 2nd dimension is lnmax*/
  int dim2;
  int dim3;
  int i, j, k;

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    dim2 = PAW->lnmax[i];
    array[i] = malloc(dim2 * sizeof(double*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      dim3 = PAW->mmax[i];
      array[i][j] = malloc(dim3 * sizeof(double));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    dim2 = PAW->lnmax[i];
    dim3 = PAW->mmax[i];
    for (j=0;j<dim2;j++) {
	  for (k=0;k<dim3;k++) {
        array[i][j][k] = 0.0;
      }
    }
  }

  /*return array*/
  return(array);
}

double*** FreeMemory_PAWvariable_threeD(double ***array, PawAtomicData * PAW, int dim1)
{
  /*frees memory form 3D double arrays for PAW variables*/
  int i, j;
  int dim2;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    dim2 = PAW->lnmax[i];
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;

  return array;
}

int** AllocateMemory_PAWmatrix_int(int **array, PawAtomicData* PAW, int dim1) 
{
  /*allocates memory for irhoij matrices for PAW*/
  int i, j;
  int dim2;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int*));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed (1D/2)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    dim2 = PAW->nrhoij[i];
    array[i] = malloc(dim2 * sizeof(int));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed(2D/2)\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    dim2 = PAW->nrhoij[i];
    for (j=0;j<dim2;j++) {
      array[i][j] = 0;
    }
  }

  /*return the array at end of function*/
  return array;
} 

double** AllocateMemory_PAWmatrix_double(double **array, PawAtomicData* PAW, int dim1) 
{
  /*allocates memory for the Dij variables in PAW*/
  int dim2;
  int i, j;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of array*/
  array = malloc(dim1 * sizeof(double*));
  /*check allocation for every dimension*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed (1D)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    dim2 = PAW->nrhoij[i];
    array[i] = malloc(dim2 * sizeof(double));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    dim2 = PAW->nrhoij[i];
    for (j=0;j<dim2;j++) {
	  array[i][j] = 0.0;
    }
  }

  return array;
}

double*** AllocateMemory_PAWthreeD_double(double ***array, PawAtomicData* PAW, AtomicVariables * ATM)
{
  /*allocates memory for the Dij_matrix*/

  /*option (1) dim2,3=lmn_size; (2) dim2,3=lnmax*/
  int i, j, k;
  int dim1, dim2, dim3;
  int typat;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  dim1 = ATM->natom;
  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double**));
  /*check if allocation was successful*/
  if (array == NULL) {
	printf("ERROR: Memory Allocation Failed\n");
	exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
	/*initialize 2nd and 3rd dimension as lmn_size for matrix*/
	typat = ATM->typat[i];
    dim2 = PAW->lmn_size[typat]; 
    dim3 = PAW->lmn_size[typat]; 

	array[i] = malloc(dim2 * sizeof(double*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(double));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
	typat = ATM->typat[i];
    dim2 = PAW->lmn_size[typat]; 
    dim3 = PAW->lmn_size[typat]; 
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
		array[i][j][k] = 0.0;
      }
    }
  }

  /*return the array at end of function*/
  return array;

}

double*** FreeMemory_PAWthreeD_double(double ***array, PawAtomicData * PAW, AtomicVariables * ATM)
{
  int i, j;
  int dim1, dim2;
  int typat;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  dim1 = ATM->natom;
  for (i=0;i<dim1;i++) {
	typat = ATM->typat[i];
	dim2 = PAW->lmn_size[typat]; 
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;

  return array;
}

double*** AllocateMemory_PAWthreeD_double_type2(double ***array, PawAtomicData* PAW, AtomicVariables * ATM)
{
  /*allocates memory for the KE_diff_matrix used in PAW*/
  int i, j, k;
  int dim1, dim2, dim3;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  dim1 = ATM->ntypat;
  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double**));
  /*check if allocation was successful*/
  if (array == NULL) {
	printf("ERROR: Memory Allocation Failed\n");
	exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
	/*initialize 2nd and 3rd dimension as lmn_size for matrix*/
    dim2 = PAW->lnmax[i]; 
    dim3 = PAW->lnmax[i]; 

	array[i] = malloc(dim2 * sizeof(double*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(double));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*zero out memory*/
  for (i=0;i<dim1;i++) {
    dim2 = PAW->lnmax[i]; 
    dim3 = PAW->lnmax[i]; 
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
		array[i][j][k] = 0.0;
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 
double*** FreeMemory_PAWthreeD_double_type2(double ***array, PawAtomicData * PAW, AtomicVariables * ATM)
{
  /*frees memory for the KE_diff_matrix in PAW*/
  int i, j;
  int dim1, dim2;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  dim1 = ATM->ntypat;
  for (i=0;i<dim1;i++) {
	dim2 = PAW->lnmax[i]; 
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;

  return array;
}



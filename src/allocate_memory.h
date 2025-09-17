#ifndef allocate_memory_H
#define allocate_memory_H  

  void Initialize_NumberGrid(NumberGrid * NG);

  void Initialize_Symmetry(Symmetry * SYM);

  void Initialize_Wavefunction(Wavefunction * WFK);

  void Initialize_BinaryGrid(BinaryGrid * BIN);

  void Initialize_AtomicVariables(AtomicVariables * ATM);

  void Initialize_VectorIndices(VectorIndices * VECT);

  void Initialize_MottJonesConditions(MottJonesConditions * MJC);

  void Initialize_TwoTheta(TwoTheta * TTH);

  void Initialize_EnergyContribution(EnergyContribution * ECON);
 
  void Initialize_PawAtomicData(PawAtomicData * PAW);

  int* AllocateMemory_oneD_int(int *array, int dim1);

  int** AllocateMemory_twoD_int(int **array, int dim1, int dim2);

  int*** AllocateMemory_threeD_int(int ***array, int dim1, int dim2, int dim3);

  int**** AllocateMemory_fourD_int(int ****array, int dim1, int dim2, int dim3, int dim4);

  double* AllocateMemory_oneD_double(double *array, int dim1);

  double** AllocateMemory_twoD_double(double **array, int dim1, int dim2);

  double*** AllocateMemory_threeD_double(double ***array, int dim1, int dim2, int dim3);

  double**** AllocateMemory_fourD_double(double ****array, int dim1, int dim2, int dim3, int dim4);

  gsl_complex*** AllocateMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2, int dim3);

  gsl_complex**** AllocateMemory_fourD_complex(gsl_complex ****array, int dim1, int dim2, int dim3, int dim4);
 
  gsl_complex***** AllocateMemory_fiveD_complex(gsl_complex *****array, int dim1, int dim2, int dim3, int dim4, int dim5);

  int* FreeMemory_oneD_int(int *array);

  int** FreeMemory_twoD_int(int **array, int dim1);

  int*** FreeMemory_threeD_int(int ***array, int dim1, int dim2);

  int**** FreeMemory_fourD_int(int ****array, int dim1, int dim2, int dim3);

  double* FreeMemory_oneD_double(double *array);

  double** FreeMemory_twoD_double(double **array, int dim1);

  double*** FreeMemory_threeD_double(double ***array, int dim1, int dim2);

  double**** FreeMemory_fourD_double(double ****array, int dim1, int dim2, int dim3);

  gsl_complex**** FreeMemory_fourD_complex(gsl_complex ****array, int dim1, int dim2, int dim3);

  gsl_complex*** FreeMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2);

  gsl_complex***** FreeMemory_fiveD_complex(gsl_complex *****array, int dim1, int dim2, int dim3, int dim4);
  
  char** AllocateMemory_twoD_char(char **array, int dim1, int dim2);

  char** FreeMemory_twoD_char(char **array, int dim1);

  void AllocateMemory_Wavefunctions(Wavefunction* WFK); 

  void FreeMemory_Wavefunctions(Wavefunction* WFK); 

  void AllocateMemory_PAWvariable(PawAtomicData* PAW, AtomicVariables * ATM);

  void FreeMemory_PAWvariable(PawAtomicData* PAW, AtomicVariables * ATM);
  
  double** AllocateMemory_PAWvariable_twoD(double** array, PawAtomicData* PAW, int dim1, int option);

  double*** AllocateMemory_PAWvariable_threeD(double*** array, PawAtomicData* PAW, int dim1);

  double*** FreeMemory_PAWvariable_threeD(double ***array, PawAtomicData * PAW, int dim1);

  double** AllocateMemory_PAWmatrix_double(double **array, PawAtomicData* PAW, int dim1);

  int** AllocateMemory_PAWmatrix_int(int **array, PawAtomicData* PAW, int dim1);

  double*** AllocateMemory_PAWthreeD_double(double ***array, PawAtomicData* PAW, AtomicVariables * ATM);

  double*** FreeMemory_PAWthreeD_double(double ***array, PawAtomicData* PAW, AtomicVariables * ATM);

  double*** AllocateMemory_PAWthreeD_double_type2(double ***array, PawAtomicData* PAW, AtomicVariables * ATM);

  double*** FreeMemory_PAWthreeD_double_type2(double ***array, PawAtomicData * PAW, AtomicVariables * ATM);

  #endif

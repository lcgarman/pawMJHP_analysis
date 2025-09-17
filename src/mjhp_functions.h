#ifndef mjhp_functions_H 
#define mjhp_functions_H 

  void find_MJregion(VectorIndices *VECT, UnitCell *UC); 

  void concatinate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP);

  void integrate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP, AtomicVariables * ATM);

  void integrate_Erange_potential(EnergyContribution * ECON, EnergyStep * ESTP, VectorIndices * VECT);
  
#endif

#ifndef extra_functions_H
#define extra_functions_H 

  void xyz2sph(double X, double Y, double Z, double * r, double * theta, double * phi);

  double**** calculate_projector_bar(double**** projector_bar, PawAtomicData * PAW, AtomicVariables * ATM, Wavefunction * WFK, UnitCell * UC) ;

#endif

#ifndef mjhp_2theta_energy_H 
#define mjhp_2theta_energy_H 

  void mjhp_2theta_local_energy(TwoTheta * TTH, NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, EnergyStep  *ESTP, EnergyContribution *ECON);

  void mjhp_2theta_nonlocal_energy(TwoTheta* TTH, PawAtomicData * PAW, AtomicVariables * ATM, UnitCell * UC, Wavefunction * WFK,EnergyStep  *ESTP, EnergyContribution *ECON);
  
#endif

#ifndef mjhp_HKL_energy_H 
#define mjhp_HKL_energy_H

  void mjhpHKL_local_energy(NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, VectorIndices *VECT, EnergyStep  *ESTP, EnergyContribution *ECON);

  void mjhpHKL_nonlocal_energy(PawAtomicData * PAW, AtomicVariables * ATM, UnitCell * UC, Wavefunction * WFK, VectorIndices *VECT, EnergyStep  *ESTP, EnergyContribution *ECON);

#endif

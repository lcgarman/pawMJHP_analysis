#ifndef boxed_band_H 
#define boxed_band_H

  typedef struct {
    double * local;
    double * nonlocal;
    double * kinetic;
    double * KE_correction;
    double * total;
    int max_Gbin;
  } BoxedBand;

  void Initialize_BoxedBand(BoxedBand * BOX); 

  void Box_Allocation(BoxedBand *BOX, NumberGrid * GRD, UnitCell * UC);

  void print_boxed_band(char filename[200], BoxedBand * BOX);

  void boxed_local(NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, BoxedBand * BOX);

  void boxed_nonlocal(PawAtomicData * PAW, AtomicVariables * ATM, UnitCell * UC, Wavefunction * WFK, BoxedBand *BOX);

  void boxed_kinetic(NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BoxedBand *BOX);

#endif

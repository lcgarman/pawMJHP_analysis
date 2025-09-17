#ifndef fftw_functions_H
#define fftw_functions_H 

  void FFTon_RealGrid(BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC);

  void FFTon_ReciprocalGrid(BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC, Symmetry* SYM);

#endif

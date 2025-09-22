#ifndef xsf_functions_H 
#define xsf_functions_H

  void print_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* BIN, AtomicVariables* ATM);

  void read_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* IND, AtomicVariables* ATM);

  void combine_xsf(NumberGrid* GRD, BinaryGrid* BIN, BinaryGrid* IND);

#endif

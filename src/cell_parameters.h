#ifndef cell_parameters_H
#define cell_parameters_H 

  void Determine_CellParameters( UnitCell* UC, NumberGrid* GRD); 

  void Calculate_FermiRadius(UnitCell* UC, FermiSphere* FS); 

  void Calculate_FermiDegree(UnitCell* UC, FermiSphere* FS); 

#endif

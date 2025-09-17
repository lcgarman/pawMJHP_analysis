#ifndef twotheta_functions_H
#define twotheta_functions_H  

  void calculate_powder_pattern(TwoTheta *TTH, BinaryGrid *BIN, NumberGrid *GRD, UnitCell *UC, Symmetry *SYM); 
  
  void fold_reflections_toBZ(TwoTheta *TTH); 
  
  void symmetry_folded_reflections(TwoTheta * TTH, Symmetry * SYM);

  void print_reflections(char filename[200], TwoTheta * TTH, FermiSphere * FS);

  void read_reflections(char filename[200], TwoTheta * TTH);

  void concatinate_twotheta_potential(EnergyContribution * ECON, EnergyStep * ESTP, TwoTheta *TTH); 

  void search_reflection(char filename[200], MottJonesConditions * MJC, double min_twotheta, double max_twotheta);

#endif

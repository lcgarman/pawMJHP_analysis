#ifndef structures_H
#define structures_H 

#include <gsl/gsl_complex.h>

typedef struct {
  char ABOfilename[100];
  char MJOUTfilename[100];
} FileCabinet;

typedef struct { 
  double voxelV;
  double bohr_cellV;
  double ang_cellV;
  double bohr_ax, bohr_ay, bohr_az;
  double bohr_bx, bohr_by, bohr_bz;
  double bohr_cx, bohr_cy, bohr_cz;
  double bohr_ax_star, bohr_ay_star, bohr_az_star; 
  double bohr_bx_star, bohr_by_star, bohr_bz_star; 
  double bohr_cx_star, bohr_cy_star, bohr_cz_star; 
  double ang_ax, ang_ay, ang_az;
  double ang_bx, ang_by, ang_bz;
  double ang_cx, ang_cy, ang_cz;
  double ang_ax_star, ang_ay_star, ang_az_star; 
  double ang_bx_star, ang_by_star, ang_bz_star; 
  double ang_cx_star, ang_cy_star, ang_cz_star;
} UnitCell; 

typedef struct { 
  int NGX, NGY, NGZ;
  int ngfftx, ngffty, ngfftz;
  int*** h_grid;
  int*** k_grid;
  int*** l_grid;
  int hmax, kmax, lmax;
} NumberGrid; 

typedef struct {
  double vec;
  double rad_Ef;
  double two_theta;
} FermiSphere; 

typedef struct {
  double bandE_min;
  double bandE_max;
  int nEstep;
  int dE_zero;
} EnergyStep;

typedef struct { 
  int nsym;
  int*** symrel;
  int mult_tot; 
  int* mult; 
  double** tnons;
} Symmetry;

typedef struct {
  double fermi;
  int nkpt;
  double* wtk;
  int nband;
  int* npw;
  double** kpt;
  int*** kg;
  double** eigen;
  double** occ;
  double**** cg;
} Wavefunction;

typedef struct {
  double*** real_grid;
  gsl_complex*** rec_grid;
  gsl_complex*** cc_rec_grid;
} BinaryGrid;

typedef struct {
  double** xred;
  int* typat;
  int ntypat;
  int natom;
  int* atomicno;
  double* Xcart;
  double* Ycart;
  double* Zcart;
} AtomicVariables;

typedef struct {
  int nHKL;
  int H;
  int K;
  int L;
  int* H_arr;
  int* K_arr;
  int* L_arr;
  int* H_posarr;
  int* K_posarr;
  int* L_posarr;
  double minr;
  double maxr;
  double scanE_start;
  double scanE_stop;
} VectorIndices;

typedef struct {
  char lattice[10];
  int ndts;
  int* jzH;
  int* jzK;
  int* jzL;
  int nH, nK, nL;
  int ngrid;
  double delta;
  int njzk;
  int nsym_jzk;
  double** JZkpt;
  double** JZk_sym;
  int* JZmult;
  double* JZwtk;
  double** BZk;
  int* Hpw;
  int* Kpw;
  int* Lpw;
  double* scanE_min;
  double* scanE_mid;
  double* scanE_max;
} MottJonesConditions;

typedef struct {
  int nrflc;
  double* two_theta;
  double* d_hkl;
  double* F_hkl;
  int* rflc_H;
  int* rflc_K;
  int* rflc_L;
  int* hpw;
  int* kpw;
  int* lpw;
  int* rflc_mult;
  int* rflc_H_sym;
  int* rflc_K_sym;
  int* rflc_L_sym;
  int* hpw_sym;
  int* kpw_sym;
  int* lpw_sym;
  double* kx;
  double* ky;
  double* kz;
  int nsym_bzk;
  double* kx_sym;
  double* ky_sym;
  double* kz_sym;
  double** BZkpt;
  double** BZkpt_sym;
} TwoTheta;

typedef struct {
  double* local;
  double* nonlocal;
  double* kinetic;
  double* KE_correction;
  double* total_potential;
  double** rflc_local;
  double** rflc_nonlocal;
  double** rflc_KE_correction;
  double** rflc_total;
} EnergyContribution;

typedef struct {
  char *search;
  char *replace;
} Modification;

typedef struct {
  int* nrhoij;
  int** irhoij;
  double** rhoij;
  double*** rhoij_matrix;
  double** Dij;
  double*** Dij_matrix;
  double***** Qhat;
  int* lnmax;
  int* lmn_size;
  int** l_orbital;
  int* mmax;
  char** psp_path;
  double* paw_rad;
  double** paw_lnrad;
  double** values;
  double** derivatives;
  double* shape_fxn_rc;
  double** ae_core_density;
  double** pseudo_core_density;
  double** pseudo_valence_density;
  double** zero_potential;
  double** local_ionic_potential;
  double*** ae_partial_wave;
  double*** pseudo_partial_wave;
  double*** projector_function;
  double** kinetic_energy_differences;
  double*** KE_diff_matrix;
  double** exact_exchange_matrix;
  double** exchange_core_core;
  gsl_complex**** c_squig;
}  PawAtomicData;

#endif

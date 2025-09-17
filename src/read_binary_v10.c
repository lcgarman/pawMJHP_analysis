#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"

void read_binary_v10(char filename[200], int option, UnitCell* UC, NumberGrid* GRD, Symmetry* SYM, Wavefunction* WFK, BinaryGrid* BIN, AtomicVariables* ATM, PawAtomicData* PAW) 
{
  /*Reads the Abinit Binary outbut files*/

  FILE * fab; /*Pointer to file fab to read Abinit output file*/
  int k, j; 
  char codvsn [8]; /*abinit version number*/
  char title [132]; /*info about psp*/
  int headform; /*header format version*/
  int fform; 
  int bandtot; /*nkpt*nband*/
  int intxc; /*grid breakdown for exchange-correlation E*/
  int date; /*date of calculation*/
  int ixc; /*index of exchange correlation function*/
  int natom; /*number of atoms*/
  int ngfftx; /*number of grid points for fft*/
  int ngffty;
  int ngfftz;
  int NGX, NGY, NGZ;
  int nkpt; /*number of grid points for kpoint generation*/
  int nspden; /*number of spin-density components*/
  int nspinor; /*number of spinorial components of the wavefunctions*/
  int nsppol; /*number of spin polatization*/
  int nsym; /*number of symmetry operations*/
  int npsp; /*number of pseudopotentials to be read*/
  int ntypat; /*number of types of atoms*/
  int pertcase; /*perturbative DFT?*/
  int usepaw; /*use projector augmented waves method*/
  double ecut; /*energy cutoff*/
  double ecutdg; /*energy cutoff for second grid in PAW*/
  double ecutsm; /*energy cutoff smearing*/
  double ecut_eff; 
  double qptnx; /*q-point re-normalized*/
  double qptny;
  double qptnz;
  double rprimd_ax; /*real space primitive translations*/
  double rprimd_ay;
  double rprimd_az;
  double rprimd_bx;
  double rprimd_by;
  double rprimd_bz;
  double rprimd_cx;
  double rprimd_cy;
  double rprimd_cz;
  double stmbias; /*scanning tunneling microscopy bias voltage*/
  double tphysel; /*temperature (physical) of electrons*/
  double tsmear; /*temp of smearing*/
  double znuclpsp; 
  double zionpsp;
  int pspso; /*spin-orbit coupling?*/
  int pspdat; /*revision date*/
  int pspcod; 
  int pspxc; /*XC fxnal for psp*/
  int lmn_size; /*spherical harmonics?*/
  int usewvl; /*use wavelet basis set*/
  int istwfkv;
  int nbandv;
  int npwv;
  int so_pspv; /*spint-orbit treatment for each psp should be by NATOM*/
  int symafmv; /*symmetries, anti-ferromagnetic characteristics by MAX_SYM*/
  int symrel_mx, symrel_my, symrel_mz; /*symmetry in real space*/
  int symrel_nx, symrel_ny, symrel_nz;
  int symrel_px, symrel_py, symrel_pz;
  int typatv; /*number of types of atoms*/
  int type; /*type of atom*/
  double kptx, kpty, kptz;
  double occopt; /*occupation numbers*/
  double tnons_x, tnons_y, tnons_z; /*translation non-symmorphic vect (MAXSYM)*/
  double znucltypatv; /*atomic number for nucs by NATOMS*/
  double* znucltypat; /*atomic number for nucs by NATOMS*/
  double wtkv;
  double residm; /*residual density*/
  double x,y,z; /*atomic coordinates*/
  double etotal; /*total E*/
  double fermi; /*fermi E as determined by the calc*/
  
  /*variables for binary real space grid*/
  int jx, jy, jz; 
  double bin_grid;

  /*variables for symmetry allocation*/
  double multiplicity;
  int mult_tot;

  /*variables for reading wavefunction*/
  int kptno; /*kpt index*/
  int npw; /*number of pw*/
  int nband; /*number of bands*/
  int pw; /*pw index*/
  int kx, ky, kz; /*reduced kpt coordinates*/
  int band; /*band index*/
  double eigenv; /*eigenvalues or band energy*/
  double occv; /*band occupancy*/
  double cgv; /*wavefunction coefficient*/

  int nshiftk_orig;
  int nshiftk; /*number of shifted grids to generate full kpt grid*/
  int mband;
  double amuv;
  int kptopt; /*kpt optimization*/
  int pawcpxocc;
  double nelect;
  double cellcharge;
  int icoulomb;
  int kptrlatt_ax; /*defines real space super-lattice corresponding to kpoint latt in rec space*/
  int kptrlatt_ay;
  int kptrlatt_az;
  int kptrlatt_bx;
  int kptrlatt_by;
  int kptrlatt_bz;
  int kptrlatt_cx;
  int kptrlatt_cy;
  int kptrlatt_cz;
  int kptrlatt_orig_ax;
  int kptrlatt_orig_ay;
  int kptrlatt_orig_az;
  int kptrlatt_orig_bx;
  int kptrlatt_orig_by;
  int kptrlatt_orig_bz;
  int kptrlatt_orig_cx;
  int kptrlatt_orig_cy;
  int kptrlatt_orig_cz;
  double shiftk_orig_x;
  double shiftk_orig_y;
  double shiftk_orig_z;
  double shiftk_x;
  double shiftk_y;
  double shiftk_z;
  char md5_pseudos[32];
  int iatom;
  int ispden;
  int cplex; /*number of complex components for potential*/
  int extra_int;
  int n;
  int nrhoijv; /*number of non-zero elements in rhoij matrix*/
  int nrhoij;
  int irhoijv; /*index of the non-zero elements in the rhoij matrix*/
  double rhoijv; /*element in rhoij matrix*/


  printf("\nReading %s\n.", filename);
  /*open abinit *_o_* output file*/
  fab=fopen(filename,"rb+"); 
  if(fab==NULL) {
    printf("ERROR: %s not found. \n", filename);
    exit(0);
  } 

  /*Initializing variables and arrays*/
  j=0; 
  mult_tot = 0;
  znucltypat = NULL;

  //BEGIN READING HEADER
  //First Block
  fread(&j, sizeof(int), 1, fab);
  j=fread(codvsn, sizeof(char), 8, fab);
  fread(&headform, sizeof(int), 1, fab);
  fread(&fform, sizeof(int), 1, fab);
  fread(&j, sizeof(int), 1, fab);

  //Second Block
  fread(&j, sizeof(int), 1, fab);
  fread(&bandtot, sizeof(int), 1, fab);
  fread(&date, sizeof(int), 1, fab);
  fread(&intxc, sizeof(int), 1, fab);
  fread(&ixc, sizeof(int), 1, fab);
  fread(&natom, sizeof(int), 1, fab);
    ATM->natom = natom;
  fread(&ngfftx, sizeof(int), 1, fab);
  fread(&ngffty, sizeof(int), 1, fab);
  fread(&ngfftz, sizeof(int), 1, fab);
  /*Storing ngfft into grid info*/
	GRD->ngfftx=ngfftx;
	GRD->ngffty=ngffty;
	GRD->ngfftz=ngfftz;
	NGX=ngfftx+1;
	NGY=ngffty+1;
	NGZ=ngfftz+1;
	GRD->NGX=NGX;
	GRD->NGY=NGY;
	GRD->NGZ=NGZ;

  fread(&nkpt, sizeof(int), 1, fab);
    WFK->nkpt = nkpt;
  fread(&nspden, sizeof(int), 1, fab);
  fread(&nspinor, sizeof(int), 1, fab);
  fread(&nsppol, sizeof(int), 1, fab);
  fread(&nsym, sizeof(int), 1, fab);
    SYM->nsym = nsym;
  fread(&npsp, sizeof(int), 1, fab);
  fread(&ntypat, sizeof(int), 1, fab);
    ATM->ntypat = ntypat;
  fread(&occopt, sizeof(int), 1, fab);
  fread(&pertcase, sizeof(int), 1, fab);
  fread(&usepaw, sizeof(int), 1, fab);
  fread(&ecut, sizeof(double), 1, fab);
  fread(&ecutdg, sizeof(double), 1, fab);
  fread(&ecutsm, sizeof(double), 1, fab);
  fread(&ecut_eff, sizeof(double), 1, fab);
  fread(&qptnx, sizeof(double), 1, fab);
  fread(&qptny, sizeof(double), 1, fab);
  fread(&qptnz, sizeof(double), 1, fab);
  fread(&rprimd_ax, sizeof(double), 1, fab);
  fread(&rprimd_ay, sizeof(double), 1, fab);
  fread(&rprimd_az, sizeof(double), 1, fab);
  fread(&rprimd_bx, sizeof(double), 1, fab);
  fread(&rprimd_by, sizeof(double), 1, fab);
  fread(&rprimd_bz, sizeof(double), 1, fab);
  fread(&rprimd_cx, sizeof(double), 1, fab);
  fread(&rprimd_cy, sizeof(double), 1, fab);
  fread(&rprimd_cz, sizeof(double), 1, fab);
  /*store primitive vectors in unitcell structure*/
	UC->bohr_ax = rprimd_ax;
	UC->bohr_bx = rprimd_bx;
	UC->bohr_cx = rprimd_cx;
	UC->bohr_ay = rprimd_ay;
	UC->bohr_by = rprimd_by;
	UC->bohr_cy = rprimd_cy;
	UC->bohr_az = rprimd_az;
	UC->bohr_bz = rprimd_bz;
	UC->bohr_cz = rprimd_cz;

  fread(&stmbias, sizeof(double), 1, fab);
  fread(&tphysel, sizeof(double), 1, fab);
  fread(&tsmear, sizeof(double), 1, fab);
  fread(&usewvl, sizeof(int), 1, fab);
  fread(&nshiftk_orig, sizeof(int), 1, fab);
  fread(&nshiftk, sizeof(int), 1, fab);
  fread(&mband, sizeof(int), 1, fab);
  fread(&j, sizeof(int), 1, fab);

  /*allocate variables for next section of read in*/
  printf( "\tAllocating Memory for Header Variables...");
  WFK->npw = AllocateMemory_oneD_int(WFK->npw, nkpt);
  SYM->symrel = AllocateMemory_threeD_int(SYM->symrel, 3, 3, nsym);
  ATM->typat = AllocateMemory_oneD_int(ATM->typat, natom);
  WFK->kpt = AllocateMemory_twoD_double(WFK->kpt, 3, nkpt);
  SYM->tnons = AllocateMemory_twoD_double(SYM->tnons, 3, nsym);
  znucltypat = AllocateMemory_oneD_double(znucltypat, ntypat);
  ATM->atomicno = AllocateMemory_oneD_int(ATM->atomicno, natom);
  WFK->wtk = AllocateMemory_oneD_double(WFK->wtk, nkpt);
  SYM->mult = AllocateMemory_oneD_int(SYM->mult, nkpt);
  ATM->xred = AllocateMemory_twoD_double(ATM->xred, 3, natom);
  printf("Done.\n");
  /*end of allocation*/

  //Third Block
  fread(&j, sizeof(int), 1, fab);
  for(j=0;j<(nkpt);j++) {
    fread(&istwfkv, sizeof(int), 1, fab);
  }
  for(j=0;j<(nkpt*nsppol);j++) {
    fread(&nbandv, sizeof(int), 1, fab);
    WFK->nband = nbandv;
  }

  for(j=0;j<(nkpt);j++) {
    fread(&npwv, sizeof(int), 1, fab);
    //printf("\t k %d\t npw = %d\n", j, npwv);
    WFK->npw[j] = npwv;
  }

  for(j=0;j<(npsp);j++) {
    fread(&so_pspv, sizeof(int), 1, fab);
  }
  for(j=0;j<(nsym);j++) {
    fread(&symafmv, sizeof(int), 1, fab);
  }
  for(j=0;j<(nsym);j++) {
    fread(&symrel_mx, sizeof(int), 1, fab);
	SYM->symrel[0][0][j] = symrel_mx;
    fread(&symrel_nx, sizeof(int), 1, fab);
	SYM->symrel[1][0][j] = symrel_nx;
    fread(&symrel_px, sizeof(int), 1, fab);
	SYM->symrel[2][0][j] = symrel_px;
    fread(&symrel_my, sizeof(int), 1, fab);
	SYM->symrel[0][1][j] = symrel_my;
    fread(&symrel_ny, sizeof(int), 1, fab);
	SYM->symrel[1][1][j] = symrel_ny;
    fread(&symrel_py, sizeof(int), 1, fab);
	SYM->symrel[2][1][j] = symrel_py;
    fread(&symrel_mz, sizeof(int), 1, fab);
	SYM->symrel[0][2][j] = symrel_mz;
    fread(&symrel_nz, sizeof(int), 1, fab);
	SYM->symrel[1][2][j] = symrel_nz;
    fread(&symrel_pz, sizeof(int), 1, fab);
	SYM->symrel[2][2][j] = symrel_pz;
    //printf(" %d\n\t%d %d %d\n\t%d %d %d\n\t%d %d %d\n\n", j, SYM->symrel[0][0][j], SYM->symrel[0][1][j], SYM->symrel[0][2][j], SYM->symrel[1][0][j], SYM->symrel[1][1][j], SYM->symrel[1][2][j], SYM->symrel[2][0][j], SYM->symrel[2][1][j], SYM->symrel[2][2][j]);
  }
  for(j=0;j<(natom);j++) {
    fread(&typatv, sizeof(int), 1, fab);
      ATM->typat[j] = typatv-1;
  }

  for(k=0;k<(nkpt);k++) {
	//reduced coordinates of kpts
    fread(&kptx, sizeof(double), 1, fab);
    fread(&kpty, sizeof(double), 1, fab);
    fread(&kptz, sizeof(double), 1, fab);
	//printf("kpt %d\t %lf %lf %lf\n", k, kptx, kpty, kptz);
    WFK->kpt[0][k] = kptx;
    WFK->kpt[1][k] = kpty;
    WFK->kpt[2][k] = kptz;
  }

  for(j=0;j<(bandtot);j++) {
    fread(&occv, sizeof(double), 1, fab);
  }

  for(j=0;j<(nsym);j++) {
    fread(&tnons_x, sizeof(double), 1, fab);
	SYM->tnons[0][j]=tnons_x;
    fread(&tnons_y, sizeof(double), 1, fab);
	SYM->tnons[1][j]=tnons_y;
    fread(&tnons_z, sizeof(double), 1, fab);
	SYM->tnons[2][j]=tnons_z;
	//printf("%d\n\t%lf %lf %lf\n", j, tnons_x, tnons_y, tnons_z);
  }
  for(j=0;j<(ntypat);j++) {
    fread(&znucltypatv, sizeof(double), 1, fab);
    znucltypat[j] = znucltypatv;
  }
  /*store atomic number by atoms*/
  for (j=0;j<natom;j++) {
    type = ATM->typat[j];
    ATM->atomicno[j] = floor(znucltypat[type]);
  }
  /*Read in the kpt weight and find multiplicity of kpts*/
  for(j=0;j<(nkpt);j++) {
    fread(&wtkv, sizeof(double), 1, fab);
    WFK->wtk[j] = wtkv;
	multiplicity = (WFK->wtk[j]/WFK->wtk[0]);
    SYM->mult[j] = ceil(multiplicity);
	mult_tot += SYM->mult[j];
	//printf("kpt: %lf %lf %lf\tmult = %d\t%d\n", WFK->kpt[0][j], WFK->kpt[1][j], WFK->kpt[2][j], SYM->mult[j], mult_tot);
  }
  SYM->mult_tot = mult_tot;
  //printf("Total nkpts = %d\n",  mult_tot);

  fread(&j, sizeof(int), 1, fab);

  //Fourth Block
  fread(&j, sizeof(int), 1, fab);
  fread(&residm, sizeof(double), 1, fab);
  for (k=0;k<natom;k++) {
    fread(&x, sizeof(double), 1, fab);
    fread(&y, sizeof(double), 1, fab);
    fread(&z, sizeof(double), 1, fab);
    ATM->xred[0][k] = x;
    ATM->xred[1][k] = y;
    ATM->xred[2][k] = z;
  }
  fread(&etotal, sizeof(double), 1, fab);
  fread(&fermi, sizeof(double), 1, fab);
  WFK->fermi = fermi;
  for (k=0;k<ntypat;k++) {
    fread(&amuv, sizeof(double), 1, fab);
  }
  fread(&j, sizeof(int), 1, fab);

  //Fifth block
  fread(&j, sizeof(int), 1, fab);  
  fread(&kptopt, sizeof(int), 1, fab);
  fread(&pawcpxocc, sizeof(int), 1, fab);
  fread(&nelect, sizeof(double), 1, fab);
  fread(&cellcharge, sizeof(double), 1, fab);
  fread(&icoulomb, sizeof(int), 1, fab);
  fread(&kptrlatt_ax, sizeof(int), 1, fab);
  fread(&kptrlatt_ay, sizeof(int), 1, fab);
  fread(&kptrlatt_az, sizeof(int), 1, fab);
  fread(&kptrlatt_bx, sizeof(int), 1, fab);
  fread(&kptrlatt_by, sizeof(int), 1, fab);
  fread(&kptrlatt_bz, sizeof(int), 1, fab);
  fread(&kptrlatt_cx, sizeof(int), 1, fab);
  fread(&kptrlatt_cy, sizeof(int), 1, fab);
  fread(&kptrlatt_cz, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_ax, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_ay, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_az, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_bx, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_by, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_bz, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_cx, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_cy, sizeof(int), 1, fab);
  fread(&kptrlatt_orig_cz, sizeof(int), 1, fab);
  fread(&shiftk_orig_x, sizeof(double), 1, fab);
  fread(&shiftk_orig_y, sizeof(double), 1, fab);
  fread(&shiftk_orig_z, sizeof(double), 1, fab);
  fread(&shiftk_x, sizeof(double), 1, fab);
  fread(&shiftk_y, sizeof(double), 1, fab);
  fread(&shiftk_z, sizeof(double), 1, fab);
  fread(&j, sizeof(int), 1, fab);

  //Sixth Block
  for(k=0;k<(npsp);k++) {
    fread(&j, sizeof(int), 1, fab); 
    fread(&title, sizeof(char), 132, fab);
    fread(&znuclpsp, sizeof(double), 1, fab);
    fread(&zionpsp, sizeof(double), 1, fab);
    fread(&pspso, sizeof(int), 1, fab);
    fread(&pspdat, sizeof(int), 1, fab);
    fread(&pspcod, sizeof(int), 1, fab);
    fread(&pspxc, sizeof(int), 1, fab);
    fread(&lmn_size, sizeof(int), 1, fab);
    fread(&md5_pseudos, sizeof(char), 32, fab);
    fread(&j, sizeof(int), 1, fab);
  }

  //Seventh Block - for PAW DATA
  PAW->nrhoij = AllocateMemory_oneD_int(PAW->nrhoij, natom);

  if (nspden != 1) {
    printf("ERROR: No code written for nspden!=1. Exiting Program.\n");
    exit(0);
  }

  if (usepaw==1) {
    fread(&j, sizeof(int), 1, fab);
    for (iatom=0;iatom<natom;iatom++) {
      for (ispden=0;ispden<nspden;ispden++) {
        fread(&nrhoijv, sizeof(int), 1, fab);
        PAW->nrhoij[iatom] = nrhoijv;
//        printf("iatom=%d\t ispden=%d\t nrhoijsel = %d\n", iatom, ispden, nrhoijv);
      }
    }
    fread(&cplex, sizeof(int), 1, fab); /*number of complex comonents == 1 or wrong*/
    fread(&nspden, sizeof(int), 1, fab); /*number of spin density*/
    /*in abinit v10 have to read an additional integer here. WHY?*/
    fread(&extra_int, sizeof(int), 1, fab); /*not sure; should be equal to 1*/
    fread(&j, sizeof(int), 1, fab);

    /*Allocate memory for the indices of and values of rhoij matrix*/
    PAW->irhoij = AllocateMemory_PAWmatrix_int(PAW->irhoij, PAW, natom);
    PAW->rhoij = AllocateMemory_PAWmatrix_double(PAW->rhoij, PAW, natom);

    fread(&j, sizeof(int), 1, fab);
    for (iatom=0;iatom<natom;iatom++) {
      for (ispden=0;ispden<nspden;ispden++) {
        nrhoij = PAW->nrhoij[iatom];
        for (n=0;n<nrhoij;n++) {
          fread(&irhoijv, sizeof(int), 1, fab);
          PAW->irhoij[iatom][n] = irhoijv;
//          printf("iatom=%d\t ispden=%d\t irho=%d \t\trhoijselect = %d\t \n", iatom, ispden, n, PAW->irhoij[iatom][n]);
        }
      }
    }

    for (iatom=0;iatom<natom;iatom++) {
      for (ispden=0;ispden<nspden;ispden++) {
        nrhoij = PAW->nrhoij[iatom];
        for (n=0;n<nrhoij;n++) {
          fread(&rhoijv, sizeof(double), 1, fab);
          PAW->rhoij[iatom][n] = rhoijv;
//          printf("iatom=%d\t ispden=%d\t irho=%d \t\trhoijp = %lf\t \n", iatom, ispden, n, PAW->rhoij[iatom][n]);
        }
      }
    }
    fread(&j, sizeof(int), 1, fab);

  }  //end of PAW data

  //END of header

  /*Free local variables*/
  znucltypat = FreeMemory_oneD_double(znucltypat);
  
  /*Option=2 => Read only header then exit*/  
  if (option == 2) {
	/*free variables that were allocated during read in*/
	printf( "\t Freeing memory from header variables...");
	WFK->npw = FreeMemory_oneD_int(WFK->npw);
	SYM->symrel = FreeMemory_threeD_int(SYM->symrel, 3, 3);
	ATM->typat = FreeMemory_oneD_int(ATM->typat);
	WFK->kpt = FreeMemory_twoD_double(WFK->kpt, 3);
	SYM->tnons = FreeMemory_twoD_double(SYM->tnons, 3);
	znucltypat = FreeMemory_oneD_double(znucltypat);
	ATM->atomicno = FreeMemory_oneD_int(ATM->atomicno);
	WFK->wtk = FreeMemory_oneD_double(WFK->wtk);
	SYM->mult = FreeMemory_oneD_int(SYM->mult);
	ATM->xred = FreeMemory_twoD_double(ATM->xred, 3);
	printf("Done.\n");
	/*end of allocation*/
	return;
  }

  /*Option=1 => Read in Binary POT/DEN info*/
  else if (option == 1) {
	printf( " \tAllocating Memory for Binary Grid...");
    BIN->real_grid = AllocateMemory_threeD_double(BIN->real_grid, NGX, NGY, NGZ);
    printf("Done.\n");
    /*Begin reading in binary grid in real space*/
    fread(&j, sizeof(int), 1, fab);
    for(jz=0;jz<NGZ;jz++) {
      for(jy=0;jy<NGY;jy++) {
        for(jx=0;jx<NGX;jx++) {
          if((jx<NGX-1)&&(jy<NGY-1)&&(jz<NGZ-1)) {
            fread(&bin_grid,sizeof(double),1,fab);
            BIN->real_grid[jx][jy][jz]=bin_grid;
          }
        }
      }
    }
    fread(&j, sizeof(int), 1, fab);
    /*END of reading in information from binary file*/
    fclose(fab);

    /*Now unwrap grid so point (jz=ngfftx)==(jz=0)*/
    for(jz=0;jz<NGZ;jz++) {
      for(jy=0;jy<NGY;jy++) {
        for(jx=0;jx<NGX;jx++) {
          if(jx==NGX-1) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][jy][jz];
          if(jy==NGY-1) BIN->real_grid[jx][jy][jz]=BIN->real_grid[jx][0][jz];
          if(jz==NGZ-1) BIN->real_grid[jx][jy][jz]=BIN->real_grid[jx][jy][0];
          if((jx==NGX-1)&&(jy==NGY-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][0][jz];
          if((jx==NGX-1)&&(jz==NGZ-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][jy][0];
          if((jy==NGY-1)&&(jz==NGZ-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[jx][0][0];
          if((jy==NGY-1)&&(jz==NGZ-1)&&(jx==NGX-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][0][0];
		}
	  }
	}
	/*End of reading in and manipulating Binary Real space real_grid*/
    printf( "Finished reading %s.\n", filename);
  } //END of if option==1

  /*if option!=1 => Read in WFK file*/
  else {

	/*Allocate Memory for wavefunction variables*/
    printf( "\tAllocating Memory for Wavefunction Variables...");
    WFK->eigen = AllocateMemory_twoD_double(WFK->eigen, nkpt, WFK->nband);
    WFK->occ = AllocateMemory_twoD_double(WFK->occ, nkpt, WFK->nband);
    fflush(stdout);
    AllocateMemory_Wavefunctions(WFK);
    printf("Done.\n");
    
	//Being reading WFK info
    for (k=0;k<nsppol;k++) {
      for (kptno=0;kptno<nkpt;kptno++) {
		fflush(stdout);
        fread(&j, sizeof(int), 1, fab);
        fread(&npw, sizeof(int), 1, fab);
        fread(&nspinor, sizeof(int), 1, fab);
        fread(&nband, sizeof(int), 1, fab);
        fread(&j, sizeof(int), 1, fab);
        fread(&j, sizeof(int), 1, fab);
        for (pw=0;pw<npw;pw++) {
          fread(&kx, sizeof(int), 1, fab);
          fread(&ky, sizeof(int), 1, fab);
          fread(&kz, sizeof(int), 1, fab);
	  	  //kxyz = planewave reduced coords
          WFK->kg[kptno][pw][0]=kx;
          WFK->kg[kptno][pw][1]=ky;
          WFK->kg[kptno][pw][2]=kz;
        }
		fread(&j, sizeof(int), 1, fab);
        fread(&j, sizeof(int), 1, fab);
        for (band=0;band<nband;band++) {
          fread(&eigenv, sizeof(double), 1, fab);
	  	  WFK->eigen[kptno][band]=eigenv;
		}
        for (band=0;band<nband;band++) {
          fread(&occv, sizeof(double), 1, fab);
		  WFK->occ[kptno][band]=occv;
        }
		fread(&j, sizeof(int), 1, fab);

		for(band=0;band<nband;band++) {
		  fread(&j, sizeof(int), 1, fab);
          for(pw=0;pw<npw;pw++) {
            fread(&cgv, sizeof(double), 1, fab);
            WFK->cg[kptno][band][pw][0]=cgv;
            fread(&cgv, sizeof(double), 1, fab);
            WFK->cg[kptno][band][pw][1]=cgv;
		  }
		  fread(&j, sizeof(int), 1, fab);

		} //END nband loop 
	  } //END nkpt loop
	}  //END nspoll loop

	/*END of reading in WFK file*/
  fclose(fab);
  printf( "Finished reading %s.\n", filename);
  } //END of else stmt for WFK files

} //END Read_About function


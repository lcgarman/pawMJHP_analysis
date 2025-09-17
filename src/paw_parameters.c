#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "paw_parameters.h"

void read_about_pseudopotential(char filename[200], PawAtomicData * PAW, AtomicVariables * ATM)
{
  /*This function reads in the abinit.out file to extract information related to the paw
 * atomicdata sets*/
  FILE * fabout; /*abinit output file string*/
  int stop; /*integer to stop the while loop*/
  int check; /*integer checking if end of file is reached*/
  char str[MAX_STRING]; /*arbitrary string*/
  char strj[MAX_STRING]; /*another arbitrary string*/
  int ntypat; /*number of types of atoms*/
  int atom_type; /*atom type*/
  char psp_path[MAX_STRING]; /*path to the paw atomicdata set*/
  int lnmax; /*number of partial waves*/
  int l_orbital; /*angular momentum quantum number*/
  int lmn_size; /*number of different l and m values*/
  int mmax; /*number of radial grid points*/
  int natom; /*number of atoms*/
  int at; /*index for atom*/
  int atomno; /*read in of atom number*/
  int typat; /*type of atom*/
  int iDij; /*index for nonlocal potential matrix*/
  double Dij_value; /*nonlocal potential matrix coefficient*/
  int i, j; /*matrix element indices*/
  int nrhoij; /*nonlocal density matrix coefficient*/
  int n; /*index for nrhoij elements*/
  int ln; /*index for lnmax values*/
  
  /*initialize variables*/
  ntypat = ATM->ntypat;
  natom = ATM->natom;

  /*allocate memory for paw parameters*/
  PAW->lnmax = AllocateMemory_oneD_int(PAW->lnmax, ntypat);
  PAW->lmn_size = AllocateMemory_oneD_int(PAW->lmn_size, ntypat);
  PAW->mmax = AllocateMemory_oneD_int(PAW->mmax, ntypat);
  PAW->psp_path = AllocateMemory_twoD_char(PAW->psp_path, ntypat, MAX_STRING);
  PAW->Dij = AllocateMemory_PAWmatrix_double(PAW->Dij, PAW, natom);

  /*manually allocate for l_orbital*/
  PAW->l_orbital = malloc(ntypat * sizeof(int*));
  if (PAW->l_orbital==NULL) {
    printf("ERROR: Memory Allocation Failed\n");
  }

  printf( "\nReading %s for pseudopotential information.\n", filename);
  /*open abinit *.out file in read mode*/
  fabout = fopen(filename, "r");
  if(fabout==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  /*read output file until end of file (EOF) is reached*/
  stop = 0;
  while (stop == 0) {
    check=fscanf(fabout, "%s", str);
    if (check==EOF) {
      stop = 1;
    }

    /*search for the word pspini in output file - begining of paw parameter info*/
    if (strcmp(str, "pspini:")==0) {
      fscanf(fabout, "%s %s %d %s %s %s %s", strj, strj, &atom_type, strj, strj, strj, psp_path);
      /*decrease atomtype value by one*/
      atom_type--;  
      /*store the path to the atomic datasets*/
      strcpy(PAW->psp_path[atom_type], psp_path);
    }

    /*search and store lnmax variable*/
    if (strcmp(str, "(lnmax)=")==0) {
      fscanf(fabout, "%d", &lnmax);
      PAW->lnmax[atom_type] = lnmax;
    }

    /*search and store lmn_size variable*/
    if (strcmp(str, "(lmn_size=")==0) {
      fscanf(fabout, "%d", &lmn_size);
      PAW->lmn_size[atom_type] = lmn_size;
    }

    /*search for orbitals then store l-values*/
    if (strcmp(str, "orbitals=")==0) {
	  /*allocate memory for second dimension of l_orbital*/ 
	  PAW->l_orbital[atom_type] = malloc(lnmax * sizeof(int));
	  /*check allocation*/
	  if (PAW->l_orbital[atom_type] == NULL) {
		printf("ERROR: Memory Allocation Failed(2D/2)\n");
		/*free previously allocated memory*/
		for (j=0;j<atom_type;j++) {
		  free(PAW->l_orbital[j]);
		}
		free(PAW->l_orbital);
		exit(0);
	  }
	  /*zero out memory*/
	  for (j=0;j<lnmax;j++) {
		PAW->l_orbital[atom_type][j] = 0;
	  }
      /*now read andstore the PAW->l_orbital values*/
	  for (ln=0;ln<lnmax;ln++) {
		fscanf(fabout, "%d", &l_orbital);
		PAW->l_orbital[atom_type][ln] = l_orbital;
	  }
    }

    /*search and store radial grid values mmax*/
    if (strcmp(str, "mmax=")==0) {
      fscanf(fabout, "%d", &mmax);
      PAW->mmax[atom_type] = mmax;
    }
    
    /*now search for Dij matrix*/
    if (strcmp(str, "Dij")==0) {
      /*allocate memory for Dij full matrix*/
      PAW->Dij_matrix =  AllocateMemory_PAWthreeD_double(PAW->Dij_matrix, PAW, ATM);

      /*read in arbitrary string*/
      fscanf(fabout, "%s", strj);
      /*look over the number of atoms, Dij and Rhoij matrix for each atom*/
      for (at=0;at<natom;at++) {
        /*initialize matrix dependent varibales*/
        iDij = 0;
        nrhoij = PAW->nrhoij[at];
		typat = ATM->typat[at];
		lmn_size = PAW->lmn_size[typat]; 
        /*if there is more than one atom, read information about atomnumber*/
        if (natom>1) fscanf(fabout, "%s %s %d", strj, strj, &atomno);

        /*matrix organized by partial wave number and m value*/
        for (i=0;i<lmn_size;i++) {
          for (j=0;j<lmn_size;j++) {
            /*read in all values of the Dij matrix*/
            fscanf(fabout, "%lf", &Dij_value);
            PAW->Dij_matrix[at][i][j] = Dij_value;
            /*only interested in half of the matrix - symmetric*/
            if (j<=i) {
              iDij++;
              /*loop over the number of rhoij elements*/
              for (n=0;n<nrhoij;n++) {
                /*if the iDij index matches any rhoij index store the value*/
                if (iDij==PAW->irhoij[at][n]) {
                  PAW->Dij[at][n] = Dij_value;
                }
              } /*end of n->nrhoij*/
			} 
		  } /*end of j->lmn_size*/
		} /* end of i->lmn_size*/
        /*read additional line at end of matrix read in*/ 
        fscanf(fabout, "%s %s %s %s %s %s", strj, strj, strj, strj, strj, strj);
	  } /*end of iatom->natom*/
    } /*end of search for Dij*/
        
  } /*end of while loop*/


  /*print out rhoij and Dij matrix elements*/
  for (at=0;at<natom;at++) {
    printf("\tatom: %d\n", at);
    nrhoij = PAW->nrhoij[at];
    for (n=0;n<nrhoij;n++) {
	  printf("\trhoij=%lf\tDij=%lf\n", PAW->rhoij[at][n], PAW->Dij[at][n]);
    }
  }

  /*close the abinit output file*/
  fclose(fabout);

} //END of  read_about_pseudopotential


void read_atomic_data(PawAtomicData * PAW, AtomicVariables * ATM)
{
  /*This function reads the atomic dataset information used in an abinit calculation*/
  FILE * fatdat; /*pointer to the atomicdata set file*/
  char psp_path[MAX_STRING]; /*path to the atomic dataset */
  int stop; /*flag integer to stop read in*/
  int check; /*integer checking for end of file*/
  int ntypat; /*number of type of atoms*/
  int typat; /*index for atom type*/
  char str[MAX_STRING]; /*arbitrary string*/
  char strj[MAX_STRING]; /*another arbitrary string*/
  char charj[3]; /*arbitrary characters*/
  char ln_label[MAX_STRING]; /*label for the partial wave*/
  int lnmax; /*number of partial waves*/
  int mmax; /*number of radial grid points*/
  int m; /*radial grid index*/
  int ln; /*partial wave index*/
  int i, j; 

  double ae_energy_kinetic; /*all-electron kinetic energy*/
  double ae_energy_xc; /*all-electron exchange and correlation energy*/
  double ae_energy_electrostatic; /*all-electron electrostatic energy*/
  double ae_energy_total; /*total all-electron energy*/
  double core_energy_kinetic; /*core kinetic energy*/
  double paw_radius; /*paw radius*/
  double ln_radius; /*radius for partial waves*/
  double occ_basis; /*occupancy of partial waves*/
  int n_quantno; /*quantum number of partial waves*/

  double values; /*radial grid values*/
  double derivatives; /*derivatives of radial grid values*/
  double shape_function_rc; /*radius for shape function*/
  double ae_core_density_rc; /*ae core density radius*/
  double ae_core_density; /*ae core density*/
  double pseudo_core_density_rc; /*ps core density radius*/
  double pseudo_core_density; /*ae core density*/
  double pseudo_valence_density_rc; /*ps valence density radius*/
  double pseudo_valence_density; /*ps valence density*/
  double zero_potential_rc; /*zero potential radius*/
  double zero_potential; /*zero potential values*/
  double local_ionic_potential_rc; /*local ionic potential radius*/
  double local_ionic_potential; /*local potential values*/
  double ae_partial_wave; /*ae partial wave values*/
  double pseudo_partial_wave; /*ps partial wave values*/
  double projector_function; /*projector function values*/
  double kinetic_energy_differences; /*kinetic energy difference values*/
  double exact_exchange_matrix; /*exact exchange values*/
  double exchange_core_core; /*exchange core values*/
  
  /*initialize variables*/
  ntypat = ATM->ntypat;
  /*allocate memory for paw atomic variables*/
  AllocateMemory_PAWvariable(PAW, ATM);

  /*read in atomic data file for each type of atom*/
  for (typat=0;typat<ntypat;typat++) {
  
    /*initialize variables for this atom type from struct*/
    strcpy(psp_path, PAW->psp_path[typat]);
    lnmax = PAW->lnmax[typat];
    mmax = PAW->mmax[typat];

    /*open atomic data file for first atom type*/
    printf( "\nReading %s for atomic data.\n", psp_path);
    fatdat = fopen(psp_path, "r");
    if(fatdat==NULL) {
      printf("%s not found. \n", psp_path);
      exit(0);
    }  

    /*set indices and flag to zero prior to read in*/
	ln=0;
    stop = 0;
    while (stop == 0) {
      check=fscanf(fatdat, "%s", str);
      if (check==EOF) {
        stop = 1;
      }

	  /*read in ae_energy*/
      if (strcmp(str, "<ae_energy")==0) {
        fscanf(fatdat, "%s %lf %c %s %lf %c %s %lf %c %s %lf %c %s", strj, &ae_energy_kinetic, charj, strj, &ae_energy_xc, charj, strj, &ae_energy_electrostatic, charj, strj, &ae_energy_total, charj, strj);
      }

	  /*read in core_energy */
      if (strcmp(str, "<core_energy")==0) {
	    fscanf(fatdat, "%s %lf %c %s", strj, &core_energy_kinetic, charj, strj);
      }

	  /*read in paw_radius*/
      if (strcmp(str, "<paw_radius")==0) {
        fscanf(fatdat, "%s %lf %s ", strj, &paw_radius, strj);
        PAW->paw_rad[typat] = paw_radius;
      }

      /*read in valence_states information*/
      if (strcmp(str, "<valence_states>")==0) {
		for (ln=0;ln<lnmax;ln++) {
          fscanf(fatdat, "%s %s %d %c %s %s %lf %c %s %lf %c %s %s %s", strj, strj, &n_quantno, charj, strj, strj, &occ_basis, charj, strj, &ln_radius, charj, strj, ln_label, strj);
          PAW->paw_lnrad[typat][ln] = ln_radius;
        }
      }

      /*read in radial grid r values*/
      if (strcmp(str, "<values>")==0) {
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &values);
          PAW->values[typat][m] = values;
		}
	  }

      /*read in derivatives of radial grid dr */
      if (strcmp(str, "<derivatives>")==0) {
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &derivatives);
          PAW->derivatives[typat][m] = derivatives;
		}
	  }

      /*reaad in sphape function*/
      if (strcmp(str, "<shape_function")==0) {
        fscanf(fatdat, "%s %s %lf %c %s", strj, strj, &shape_function_rc, charj, strj);
        PAW->shape_fxn_rc[typat] = shape_function_rc;
      }

      /*read in all electron core density*/
      if (strcmp(str, "<ae_core_density")==0) {
        fscanf(fatdat, "%s %s %lf %c %s", strj, strj, &ae_core_density_rc, charj, strj);
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &ae_core_density);
          PAW->ae_core_density[typat][m] = ae_core_density;
        }
      }

      /*read in pseudo_core_density*/
      if (strcmp(str, "<pseudo_core_density")==0) {
        fscanf(fatdat, "%s %s %lf %c %s", strj, strj, &pseudo_core_density_rc, charj, strj);
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &pseudo_core_density);
          PAW->pseudo_core_density[typat][m] = pseudo_core_density;
        }
      }
        
      /*read in pseudo_valence_density*/
      if (strcmp(str, "<pseudo_valence_density")==0) {
        fscanf(fatdat, "%s %s %lf %c %s", strj, strj, &pseudo_valence_density_rc, charj, strj);
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &pseudo_valence_density);
          PAW->pseudo_valence_density[typat][m] = pseudo_valence_density;
        }
      }

      /*read in  zero_potential*/
      if (strcmp(str, "<zero_potential")==0) {
        fscanf(fatdat, "%s %s %lf %c %s", strj, strj, &zero_potential_rc, charj, strj);
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &zero_potential);
          PAW->zero_potential[typat][m] = zero_potential;
        }
      }
    
      /*read in  local_ionic_potential*/
      if (strcmp(str, "<blochl_local_ionic_potential")==0) {
        fscanf(fatdat, "%s %s %lf %c %s", strj, strj, &local_ionic_potential_rc, charj, strj);
        for (m=0;m<mmax;m++) {
          fscanf(fatdat, "%lf", &local_ionic_potential);
          PAW->local_ionic_potential[typat][m] = local_ionic_potential;
        }
      }
    
      /*read in data dependend on basis number*/
	  if (strcmp(str, "</blochl_local_ionic_potential>")==0) {
        /*loop over every partial wave*/
        for (ln=0;ln<lnmax;ln++) {
		/*read in  ae_partial_wave*/
		  fscanf(fatdat, "%s %s %s %s ", strj, strj, strj, strj);
		  for (m=0;m<mmax;m++) {
			fscanf(fatdat, "%lf", &ae_partial_wave);
			PAW->ae_partial_wave[typat][ln][m] = ae_partial_wave;
		  }
          fscanf(fatdat, "%s", strj);

  		  /*read in  ps_partial_wave*/
		  fscanf(fatdat, "%s %s %s %s ", strj, strj, strj, strj);
		  for (m=0;m<mmax;m++) {
			fscanf(fatdat, "%lf", &pseudo_partial_wave);
			PAW->pseudo_partial_wave[typat][ln][m] = pseudo_partial_wave;
		  }
          fscanf(fatdat, "%s", strj);

		  /*read in  projector_functions*/
		  fscanf(fatdat, "%s %s %s %s", strj, strj, strj, strj);
		  for (m=0;m<mmax;m++) {
			fscanf(fatdat, "%lf", &projector_function);
			PAW->projector_function[typat][ln][m] = projector_function;
		  }
          fscanf(fatdat, "%s", strj);

		} /*end of loop over ln->lnmax*/
	  }  

	  /*read in  kinetic_energy_differences*/
	  if (strcmp(str, "<kinetic_energy_differences>")==0) {
        for (i=0;i<lnmax;i++) {
          for (j=0;j<lnmax;j++) {
		    fscanf(fatdat, "%lf", &kinetic_energy_differences);
		    PAW->KE_diff_matrix[typat][i][j] = kinetic_energy_differences;
          }
		}
	  }

	  /*read in  exact_exchange_matrix*/
	  if (strcmp(str, "<exact_exchange_X_matrix>")==0) {
		for (ln=0;ln<(lnmax*lnmax);ln++) {
		  fscanf(fatdat, "%lf", &exact_exchange_matrix);
		  PAW->exact_exchange_matrix[typat][ln] = exact_exchange_matrix;
		}
	  }

	  /*read in exact_exchange*/
	  if (strcmp(str, "<exact_exchange")==0) {
		fscanf(fatdat, "%s %lf %c %s", strj, &exchange_core_core, charj, strj);
	  }

	} /*end of while loop*/

    /*close this paw data set file*/
    fclose(fatdat);

    printf("\tatomic data successfully read in for atom type %d.\n", typat);
  } /*end of loop over ntypat*/
  
} //read_atomic_data


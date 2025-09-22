#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "xsf_functions.h"

void print_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* BIN, AtomicVariables* ATM)
{
  /*This function prints the density information in an XSF file format.*/
  FILE * fxsf; /*pointer to the xsf file*/
  int jx, jy, jz; /*indices for direct space number grid*/
  int j; /*index*/
  int natom; /*number of atoms*/
  double x, y, z; /*atomic coordinates*/
  double* Xcart; /*cartesian coordinates of atoms*/
  double* Ycart;
  double* Zcart;
  int NGX, NGY, NGZ; /*direct space number grid*/
  int line_counter; /*counts how many items have been printed on a line*/

  /*initialize variables*/
  natom = ATM->natom;
  NGX = GRD->NGX;
  NGY = GRD->NGY;
  NGZ = GRD->NGZ;
  Xcart = NULL;
  Ycart = NULL;
  Zcart = NULL;

  /*open the xsf file in write mode*/
  fxsf=fopen(filename, "w");
  if (fxsf==NULL) {
    printf("ERROR: %s not found\n", filename);
    exit(0);
  }

  /*allocate mem for XYZcart*/
  Xcart = AllocateMemory_oneD_double(Xcart, natom);
  Ycart = AllocateMemory_oneD_double(Ycart, natom);
  Zcart = AllocateMemory_oneD_double(Zcart, natom);

  /*if Xcart variable already exists then fill the local arrays*/
  if (ATM->xred == NULL) {
    for (j=0;j<natom;j++) {
      Xcart[j] = ATM->Xcart[j];
      Ycart[j] = ATM->Ycart[j];
      Zcart[j] = ATM->Zcart[j];
    }
  } 
  else {
    /*if Xcart does not exist, then calculate it here*/
	for (j=0;j<natom;j++) {
	  x = ATM->xred[0][j];
	  y = ATM->xred[1][j];
	  z = ATM->xred[2][j];
	  Xcart[j] = x*UC->ang_ax+y*UC->ang_bx+z*UC->ang_cx; 
	  Ycart[j] = x*UC->ang_ay+y*UC->ang_by+z*UC->ang_cy; 
	  Zcart[j] = x*UC->ang_az+y*UC->ang_bz+z*UC->ang_cz; 
	}
  }

  printf("\nPrinting to XSF file: %s.\n", filename);

  /*begin printing header to xsf file*/
  fprintf(fxsf, " DIM-GROUP\n"); 
  fprintf(fxsf, " 3 1 \n");
  fprintf(fxsf, "PRIMVEC\n"); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_ax, UC->ang_ay, UC->ang_az); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_bx, UC->ang_by, UC->ang_bz); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_cx, UC->ang_cy, UC->ang_cz); 

  /* coordinates of primitive lattice */
  fprintf(fxsf, "PRIMCOORD\n" ); 
  fprintf(fxsf, "\t\t%d%3d\n", natom, 1); 
  for(j=0;j<natom;j++) {
    fprintf(fxsf, "%9d%20.10lf%20.10lf%20.10lf\n", ATM->atomicno[j], Xcart[j], Ycart[j], Zcart[j]);
  }
  /* coordinates of primitive lattice */
  fprintf(fxsf, "ATOMS\n" ); 
  for(j=0;j<natom;j++) {
    fprintf(fxsf, "%9d%20.10lf%20.10lf%20.10lf\n", ATM->atomicno[j], Xcart[j], Ycart[j], Zcart[j]);
  }

  /*print header before number grid block*/
  fprintf(fxsf," BEGIN_BLOCK_DATAGRID3D\n");
  fprintf(fxsf," Written_by_print_XSF\n");
  fprintf(fxsf," DATAGRID_3D_DENSITY\n");
  fprintf(fxsf, "\t\t%d\t%d\t%d\n", NGX, NGY, NGZ);
  /*shift grid is set to 0 0 0 */
  fprintf(fxsf, "%lf\t%lf\t%lf\n", 0.0, 0.0, 0.0); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_ax, UC->ang_ay, UC->ang_az); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_bx, UC->ang_by, UC->ang_bz); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_cx, UC->ang_cy, UC->ang_cz); 

  /*now loop over direct space grid and print density*/
  line_counter=0; 
  for(jz=0;jz<NGZ;jz++) {
    for(jy=0;jy<NGY;jy++) {
      for(jx=0;jx<NGX;jx++) {
        line_counter++;
        /*print density on grid*/
        fprintf(fxsf, "\t%.10lf" , BIN->real_grid[jx][jy][jz]);
        /*when line counter reaches 6 go to next line*/
        if(line_counter==6) {
          fprintf(fxsf, "\n");
          line_counter=0;
        }
      }
    }
  }
  fprintf(fxsf, "\nEND_DATAGRID_3D\n" );
  fprintf(fxsf, "END_BLOCK_DATAGRID3D" );
  
  /*close xsf file*/
  fclose(fxsf);

  printf("Done printing to %s.\n", filename);

  /*Check total density*/
  double coeff_total;
  coeff_total = 0.0;
  for(jx=0;jx<(NGX-1);jx++) {
    for(jy=0;jy<(NGY-1);jy++) {
	  for(jz=0;jz<(NGZ-1);jz++) {
        coeff_total += BIN->real_grid[jx][jy][jz];
      }
    }
  }
  //printf("\tTotal Normalized Density = %lf\n", coeff_total*UC->voxelV);
  printf("\tTotal Normalized Density = %lf\n", coeff_total);
  
  /*free allocated vars*/
  Xcart = FreeMemory_oneD_double(Xcart);
  Ycart = FreeMemory_oneD_double(Ycart);
  Zcart = FreeMemory_oneD_double(Zcart);
  BIN->real_grid = FreeMemory_threeD_double(BIN->real_grid, NGX, NGY);
}  /*END of print_XSF function*/

void read_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* IND, AtomicVariables* ATM)
{
  /*This function reads in an XSF file and stores the data on the ngfft grid*/
  FILE * fxsf;
  int jx, jy, jz; /*indices for direct space number grid*/
  int j; 
  char strv[20]; /*arbitrary string*/
  int intv; /*arbitrary integer*/
  double doubv; /*arbitrary double*/
  int NGX, NGY, NGZ; /*number grid in real space*/
  int natom; /*number of atoms*/
  int atomicnov; /*atomic number of the atoms*/
  double Xcartv, Ycartv, Zcartv; /*cartesian coordinates of the atoms*/
  double gridv; /*varibale to read grid info to*/
  double total_xsf; /*total summed value on xsf grid*/

  /*open xsf file in read mode*/
  printf("\nReading xsf file: %s.\n", filename);
  fxsf=fopen(filename, "r");
  if (fxsf==NULL) {
    printf("ERROR: %s not found\n", filename);
    exit(0);
  }
  total_xsf = 0.0;

  /*begin reading header to xsf file*/
  fscanf(fxsf, "%s", strv); 
  fscanf(fxsf, "%d %d", &intv, &intv);
  fscanf(fxsf, "%s", strv); 
  /*read in primitive lattice vectors*/
  fscanf(fxsf, "%lf %lf %lf", &UC->ang_ax, &UC->ang_ay, &UC->ang_az); 
  fscanf(fxsf, "%lf %lf %lf", &UC->ang_bx, &UC->ang_by, &UC->ang_bz); 
  fscanf(fxsf, "%lf %lf %lf", &UC->ang_cx, &UC->ang_cy, &UC->ang_cz); 
  /* coordinates of atoms in lattice */
  fscanf(fxsf, "%s", strv); 
  fscanf(fxsf, "%d %d", &natom, &intv); 
  /*store the number of atoms and allocate memory for the lattice vectors*/
  ATM->natom = natom;
  ATM->atomicno = AllocateMemory_oneD_int(ATM->atomicno, natom);
  ATM->Xcart = AllocateMemory_oneD_double(ATM->Xcart, natom);
  ATM->Ycart = AllocateMemory_oneD_double(ATM->Ycart, natom);
  ATM->Zcart = AllocateMemory_oneD_double(ATM->Zcart, natom);
  /*read in the atomic number followed by the coordinates for each atom*/
  for(j=0;j<natom;j++) {
    fscanf(fxsf, "%d %lf %lf %lf\n", &atomicnov, &Xcartv, &Ycartv, &Zcartv);
    /*store these values in array*/
    ATM->atomicno[j] = atomicnov;
    ATM->Xcart[j] = Xcartv;
    ATM->Ycart[j] = Ycartv;
    ATM->Zcart[j] = Zcartv;
  }
  
  /*read the atomic reduced coordinates in again*/
  fscanf(fxsf, "%s", strv); 
  for(j=0;j<natom;j++) {
    fscanf(fxsf, "%d %lf %lf %lf\n", &atomicnov, &Xcartv, &Ycartv, &Zcartv);
  }
  /*read in three arbitrary strings that serve as header to data grid*/
  fscanf(fxsf," %s", strv);
  fscanf(fxsf," %s", strv);
  fscanf(fxsf," %s", strv);

  /*read in the size of the numeber grid*/
  fscanf(fxsf, "\t%d\t%d\t%d\n", &NGX, &NGY, &NGZ);
  GRD->NGX = NGX;
  GRD->NGY = NGY;
  GRD->NGZ = NGZ;
  /*allocate memory for data on ngfft grid*/
  IND->real_grid = AllocateMemory_threeD_double(IND->real_grid, NGX, NGY, NGZ);
  
  /*shift grid is set to 0 0 0 */
  fscanf(fxsf, "%lf\t%lf\t%lf\n", &doubv, &doubv, &doubv ); 
  /*read in primitive lattice vectors again*/
  fscanf(fxsf, "%lf %lf %lf", &UC->ang_ax, &UC->ang_ay, &UC->ang_az); 
  fscanf(fxsf, "%lf %lf %lf", &UC->ang_bx, &UC->ang_by, &UC->ang_bz); 
  fscanf(fxsf, "%lf %lf %lf", &UC->ang_cx, &UC->ang_cy, &UC->ang_cz); 

  /*begin reading in data on ngfft grid*/
  total_xsf = 0.0;
  for(jz=0;jz<NGZ;jz++) {
    for(jy=0;jy<NGY;jy++) {
      for(jx=0;jx<NGX;jx++) {
        fscanf(fxsf, "%lf" , &gridv);
        /*store grid in a the real grid array*/
        IND->real_grid[jx][jy][jz] = gridv;
        total_xsf += gridv;
      }
    }
  }
  /*read end of file lines*/
  fscanf(fxsf, "%s", strv );
  fscanf(fxsf, "%s", strv );

  /*close the file*/
  fclose(fxsf);
 
  printf("\tTotal summed values on XSF grid = %e\n", total_xsf);
  
} 
/*END of outputXSF function*/

void combine_xsf(NumberGrid* GRD, BinaryGrid* BIN, BinaryGrid* IND)
{
  /*this function combines two xsf grids onto one grid*/
  int jx, jy, jz; /*indices for real space ngfft grid*/
  int NGX, NGY, NGZ; /*number grid in real space*/
  double  total_combo; /*total summed combined grid values*/

  /*initialize variables*/
  NGX = GRD->NGX;
  NGY = GRD->NGY;
  NGZ = GRD->NGZ;
  total_combo = 0.0;

  printf("\nConcatinating XSF grids.\n");
  /*loop over grid to concatinate the grids*/
  for (jz=0;jz<NGZ;jz++) {
	for (jy=0;jy<NGY;jy++) {
	  for (jx=0;jx<NGX;jx++) {
        BIN->real_grid[jx][jy][jz] += IND->real_grid[jx][jy][jz];
        total_combo += BIN->real_grid[jx][jy][jz];
	  }
	}
  }
        
  printf("\tTotal concatinated XSF grid = %lf\n", total_combo);

  /*free memory from the individual grid*/
  IND->real_grid = AllocateMemory_threeD_double(IND->real_grid, NGX, NGY, NGZ);
 
}

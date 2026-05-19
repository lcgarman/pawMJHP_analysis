# pawMJHP_analysis
The Projector Augmented Wave Mott-Jones Hamilton Population analysis and Mott-Jones Electron Density method. A program for analyzing the Mott-Jones effect in intermetallic compounds and for exploring how local atomic configurations can support of template the Mott-Jones effects.

## Installation
Provided you have an SSH key to clone the repository, run these commands:
```
git clone git@github.com:lcgarman/pawMJHP_analysis.git
cd pawMJHP_analysis/src
make 
```
This will clone the repo into your current directory, change to the necessary directory, and run the makefile to compile the necessary executables. You may copy the binaries into any directory in your $PATH (/usr/local/bin, for instance). 

## Dependencies

### Compilation

The GNU Scientific Library, libxc, and fftw3 are required for compiling the pawMJHP_analysis package.

### Running Calculation

This method uses the Projector Augmented Wave (PAW) formalism. PAW atomic datasets should be in the XML format. Atomic datasets for 10 elements commonly used within the Mott-Jones model are included in this repository under the MJHP_optimized_AtomicDataset directory. The MJHP_analysis expects output from abinit version-10 and the necessary binary abinit output files should be named: FILENAME_o_XXX, where XXX refers to the different file types specified below. 

I recommend working in three separate directories for the pawMJHP_analysis. Let’s call these directories:

  2theta
  
  HKL

  HKL_density

We will begin with the two theta MJHP analysis which provides a broad overview of a system's potential for Mott-Jones effects. 
Copy from the single-point energy calculation into the 2theta directory: the abinit FILENAME.in, FILENAME.files, and FILENAME_i_DEN (the output density from the single-point energy to be used as the starting ground-state density for this non-selfconsistent calculation).

One additional input file is required in this directory: FILENAME.mjin. Start by making this file in the 2theta directory with the following format:
```
IN_FILENAME
OUT_COMPOUND
vec #
xY
```
The first line, IN_FILENAME, should be the file string corresponding to your abinit binary output files. The OUT_COMPOUND is what you would like the output files from the MJHP_analysis program to be named. The third line expects a string followed by a number (double or integer) that corresponds to the number of valence electrons in the primitive unit cell of the compound. The fourth line, xY, should be two letters where the first (x) corresponds to the crystal family and the second to the centering (Y). This line is case sensitive, so x should be lowercase and Y uppercase, identical to the letters of a Pearson symbol. 

Once the FILENAME.mjin file and the necessary abinit files are within the 2theta directory run the executable:
```
prep_paw2theta FILENAME.mjin
```
This program will find the high symmetry k-points of interest and copy them into your abinit FILENAME.in file, as well as set up an OUT_COMPOUND.rflc file that is necessary for the paw_mjhp2theta analysis that will follow. Run a non-self consistent abinit calculation using the modified FILENAME.in input file. Following the completion of this job, in the same directory, run the following executable to perform the first pawMJHP_analysis:
```
paw_mjhp2theta  FILENAME.mjin >& FILENAME_mjhp2theta.log
```
The ">& FILENAME_mjhp2theta.log" is optional and will print the command line output to that log file. This executable will generate the OUT_FILENAME_2theta.mjout file that can be plotted using the matlab program:
```
plot_MJHP_2theta('OUT_COMPOUND_2theta.mjout')
```
This completes the first MJHP analysis. 

Now, lets focus on the MJHP analysis for specific families of reciprocal lattice vectors.

Copy from the single-point energy calculation into the HKL directory: the abinit FILENAME.out, FILENAME_o_POT, and FILENAME_o_WFK files. Also, copy over the FILENAME.mjin file from the 2theta directory. This file now requires a few additional lines:
```
IN_FILENAME
OUT_COMPOUND
vec #
xY
nHKL n
HKL
h1 k1 l1
h2 k2 l2
hn kn ln
```
The first four lines stay the same from the 2theta to HKL calculations. The fifth line, should be a string (arbitrary) followed by an integer of how many datasets you wish to run (n, ie 3). The sixth line should be a string “HKL”. Following this line there should be a list of the conventional HKL indices of interest. The HKL indices should be separated by a space (H K L), with each set of HKL indices placed on a new line. The number of lines following HKL on the sixth line should correspond to the integer (n) following the "nHKL" on the fifth line. Following this set up, run the executable:
```
paw_mjhpHKL FILENAME.mjin >& FILENAME_mjhpHKL.log
```
Again, the ">& FILENAME_mjhpHKL.log" is optional and will print the command line output to the log file. This executable will generate one OUT_COMPOUND_HnKnLn.mjout file for each HKL dataset. The results from this analysis may be plotted using the matlab program:
```
plot_MJHP_HKL('OUT_COMPOUND_h1k1l1.mjout', minY, maxY)
```
The  minY and maxY values correspond to the min and max energy range in eV you wish to view. 

This completes the MJHP_analysis. 

Now, lets extract the density associated with the Mott-Jones effects from these different families of reciprocal lattice vectors. 
Copy from the single-point energy calculation into the density_HKL directory: the abinit FILENAME_o_WFK file. Also, copy over the FILENAME.mjin file from the HKL directory. We again require a few additional lines:
```
IN_FILENAME
OUT_COMPOUND
vec #
xY
nHKL n
HKL
h1 k1 l1
h2 k2 l2
hn kn ln
scanE
minE1 0.0 maxE1
minE2 0.0 maxE2
minE3 0.0 maxE3
```
The lines from the HKL calculation stay the same, however, we have appended  additional lines to the bottom of this file. Following our list of HKL indices, we have a string "scanE" where the string itself is arbitrary. Following this line, there should be a list of the minimum, middle, and maximum energies to sample for the Mott-Jones effects for each HKL indice, where the first energy range corresponds to the first set of HKL indices, the second to the second HKL indices, and so on. Therefore, there should be an equal number of lines following the "HKL" string and the "scanE" string which should both match the integer "n" following the string "nHKL" on the fifth line. The Mott-Jones electron density will be extracted from the minE1 to 0.0 and then from 0.0 to maxE1 for HKL incides h1, k1, and l1. Following this set up, run the executable:
```
paw_mjhpHKL_density FILENAME.mjin >& FILENAME_mjhpDEN.log
```
Again, the ">& FILENAME_mjhpHKL.log" is optional and will print the command line output to the log file. This executable will generate two files for each HKL dataset: OUT_COMPOUND_HKL_minE_0.0.xsf and OUT_COMPOUND_HKL_0.0_maxE.xsf. The results from this analysis may be visualized using any software with XSF viewing capabilities (ie VESTA).

## Example .mjin File
Below is a sample mjin file for Cu5Zn8:
```
Cu5Zn8

Cu5Zn8_gammabrass

vec 42.0

cI

nHKL 3

HKL

3 3 0

4 1 1 

-4 1 1

scanE

-2.3 0.0 2.1

-1.9 0.0 1.6

-1.9 0.0 1.6
```

## Citation 
If you use the MJHP_analysis in your research, please cite the following papers:

(1) Garman, L. C.; Fredrickson, D. C. Interactions on the Fermi Sphere in an AllElectron World: The Projector Augmented Wave Implementation of the Mott-Jones Hamilton Population and the Three Archetypical Sphere Packings. The Journal of Physical Chemistry C 2026, under review.
(2) Garman, L. C.; Fredrickson, D. C. The Mott-Jones Electron Crystal: Patterning Atomic Positions for Pseudogap Formation in HumeRothery Phases. Inorganic Chemistry 2026, under review.


# TAURUS_vap 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4381279.svg)](https://doi.org/10.5281/zenodo.4381279)
[![License: GPL v3](https://img.shields.io/github/license/project-taurus/taurus_vap)](https://www.gnu.org/licenses/gpl-3.0)

## Presentation
We present the numerical code TAURUS_vap that solves the variation after particle-number projection (VAPNP) equations for real general Bogoliubov quasi-particle states represented in a
Spherical Harmonic Oscillator (SHO) basis. The model space considered is invariant under spatial and isospin rotations but no specific SHO basis is assumed such that the code can carry out
either valence-space or no-core calculations (depending on the hamiltonian file used as input). 
The variational procedure can be performed under several constraints on the expectation value of a variety of operators such as the multipole deformation, the pairing field or the components of the angular momentum.
In addition, no number parity is assumed for the Bogoliubov quasi-particle states such that the code can be used
to describe even-even, odd-even and odd-odd nuclei. Finally, the code can also be used to carry out simpler Hartree-Fock (HF) and Hartree-Fock-Bogoliubov (HFB) calculations.

This code is part of the larger project TAURUS that aims at developing a numerical suite centered around the concept of symmetry-projected variational calculations and is supported by the the European Union’s Horizon 2020 research and innovation programme.

## Compilation
We provide scripts to compile the code in two alternative ways. The first one is based on `Bash` scripting whereas the second one is relying on the `Make` compilation tool.
Pick the method you prefer but be aware that we cannot guarantee that the scripts, as is, will work on your system. In particular, if it does not work we encourage you to
check your compiler version (there are a few Fortran 2003/2008 commands that might not be implemented in old compilers) and the path to link the BLAS/LAPACK libraries (**required**).

### Bash script
To compile the code, go to the main directory and enter the command
```
bash compile.sh FC TH
```
where `FC` and `TH` are the two following arguments
* `FC` = `gfortran`, `ifort`, `mpiifort` or `mpif90`  
Fortran compiler.  
* `TH` = `omp` or `none`  
Option to enable OpenMP threading.

For example, to compile using the intel compiler and OpenMP, execute the command
```
bash compile.sh ifort omp
```
The script can also be launched without any argument
```
bash compile.sh
```
In such case, the script will assume the values `gfortran` and `none` as first and second argument, respectively.

If the compilation is successful, i.e. the executable file `taurus_vap.exe` was created and moved to the `exe` directory, the script will print "*compilation successful.*" as last message. 
Otherwise, the script will print the message "*compilation failed.*".
The cleaning of files and directories is automatically perfomed by the script after the compilation.

### Makefile
To compile the code, go to the main directory and enter the command
```
make all FC= TH=
```
where `FC` and `TH` are the same arguments as above. The arguments have to be entered *with* the name of the variable.

For example, to compile using the intel compiler and OpenMP, execute the command
```
make all FC=ifort TH=omp
```
The script can also be launched without any argument
```
make
```
In such case, the script will assume the recipe `all` and the values `FC=gfortran` and `TH=none`.  

To clean the directories before or after a fresh compilation, type
```
make clean
```
It will automatically remove the `*.mod` and `*.o` files in the the `mod/` and `obj/` subdirectories, respectively. 
You can also use the command `make deepclean` if you want to delete the directories `mod/` and `obj/` themselves.

## Execution
To execute the code, go in the directory containing the executable file and type the command
```
./taurus_vap.exe < input.txt
```
where `input.txt` is the STDIN file containing the inputs parameters. The details concerning the format of STDIN can be found in the file `extras/manual_input.pdf`.
The code also requires, in the same directory, the file defining the Hamiltonian (and model space), the name of which is written in the STDIN, and the various files containing its matrix elements.
The details concerning the format of Hamiltonian files can be found in the file `extras/manual_hamiltonian.pdf`.

To simplify the execution of the code, we provide the script `launch.sh` that performs all the necessary steps to run a calculation. 
To use it, go in the main directory and type the command

```
bash launch.sh 
```

At the end of a successful run, the code will print in the STDOUT the following sentence
```
This is the end, my only friend, the end.
```

### Input files
As explained above, the code will require as input files, the STDIN and the various files defining the Hamiltonian an model space. 
In addition, it is possible to provide a file containing a wave function that will be used as initial wave function in the calculation.

See the file `extras/manual_input.pdf` for more details.

### Output files
During its execution, the code prints various information (e.g. inputs and model space used, expectation value of the energy at each iteration, etc.) in the STDOUT. 
We recommend to store the printing in a file `output.txt` by typing
```
./taurus_vap.exe < input.txt > output.txt
```
or

```
bash launch.sh > ouput.txt
```

Additionally, the code will produce other files containing relevant informations such as the occupation numbers or the eigenvalues of the single-particle hamiltonian.
More importantly, the code will write the final wave function obtained at the end of the iterative procedure in a file.
The names of all the files produced during a run are recalled in the STDOUT.

See the file `extras/manual_input.pdf` for more details.

## Examples
We provide 2 examples of simple calculations in the `examples` directory:
* Calculation of the HFB/PNVAP minima for the nuclei Mg24, Mg25 and Al26 in the sd-shell using the USDB interation.
* Analysis of the rate of convergence for the gradient+momentum algorithm for Mg24 in the sd-shell with the USDB interaction.

To run an example calculation, go the appropriate subdirectory and execute the script
```
bash launch_example.sh
```
Note that the examples have to be launched using a version of the code compiled *without* MPI.
The first example takes about 3 minutes to be completed on a recent desktop computer. The second example is slightly more intensive and may take 15 to 20 minutes
before completion.

These examples are somewhat trivial but can be used to become familiar with the code and check that it runs properly.
In particular, the results we obtained performing these calculations can be found in the article associated with the publication of the code (see section Citation).

## Citation
If you use this code in your research work and publications, please cite us:

> B. Bally, Adrían Sánchez-Fernández and Tomás R. Rodríguez  
> Symmetry-projected variational calculations with the numerical suite TAURUS  
> I. Variation after particle-number projection  
> Eur. Phys. J. A 57, 69 (2021)   
> https://doi.org/10.1140/epja/s10050-021-00369-z

## Additional informations

### License
TAURUS_vap is licensed under GNU General Public License version 3 (see LICENSE.txt).

### Funding
The project TAURUS was initially supported by the the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No. 839847.  
https://cordis.europa.eu/project/id/839847

### Contributors 
For the time being, the people that contributed to the code are:
* Benjamin Bally (CEA Paris-Saclay)
* Tomás R. Rodríguez Frutos (Universidad Autónoma de Madrid)
* Adrián Sánchez-Fernández (Universidad Autónoma de Madrid)

We also thank the people that helped us benchmark and test the code in its early versions.

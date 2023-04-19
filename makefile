#################
#  Description  #
#################
 
# This is an example of makefile to compile TAURUS_vap. The code requires the
# BLAS/LAPACK libraries. When using the intel compiler "ifort", we recommend
# to use of their specific Math Kernel Library (MKL).

# This script is only given as an example and we do not guarantee that it will
# work on your system. In particular, check the version of your compiler and
# the directories of the libraries.

#################
#  Directories  #
#################

SRCDIR=src
OBJNAM=obj
OBJDIR=$(SRCDIR)/${OBJNAM}
MODDIR=$(SRCDIR)/mod
EXEDIR=exe

#############################################
#  Fortran compiler, options and libraries  #
#############################################

# default 
FC=gfortran
TH=none

ifeq ($(FC),$(filter $(FC), gfortran mpif90))
 OPT=-O3 -J$(MODDIR)
 LIB=-L/usr/lib -llapack -lblas 
 ifeq ($(TH),omp)
   OPT=-O3 -J$(MODDIR) -fopenmp
 endif
else ifeq ($(FC),$(filter $(FC), ifort mpiifort))
 OPT=-O3 -module $(MODDIR) -qmkl -lmkl_sequential
 LIB=
 ifeq ($(TH),omp)
   OPT=-O3 -module $(MODDIR) -qmkl -qopenmp -lmkl_intel_thread
 endif
endif

#######################
#  Files definitions  #
#######################

code=taurus_vap
exec=$(code).exe

wSRC77=$(wildcard $(SRCDIR)/*.f)
wSRCTMP=$(wildcard $(SRCDIR)/*.f90)

ifeq ($(FC),$(filter $(FC), ifort gfortran))
  wSRC90=$(shell echo $(wSRCTMP) | sed "s/$(SRCDIR)\/module_parallelization.f90//g")
else
  wSRC90=$(wSRCTMP)
endif

SRC77=$(shell echo $(wSRC77) | sed "s/$(SRCDIR)/$(SRCDIR)\/${OBJNAM}/g")
SRC90=$(shell echo $(wSRC90) | sed "s/$(SRCDIR)/$(SRCDIR)\/${OBJNAM}/g")

OBJ90=$(SRC90:.f90=.o)
OBJ77=$(SRC77:.f=.o)

################
#  Make rules  #
################

.PHONY: all clean deepclean

# Main file
#==========

all: $(code)

$(code): $(OBJ90) $(OBJ77) | $(OBJDIR)/ $(MODDIR)/ $(EXEDIR)/
	$(FC) $(OPT) -o $(EXEDIR)/$(exec) $^ $(LIB)
	@echo "Compilation finished."

# General rules 
#==============

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)/ $(MODDIR)/
	@cp $< $(OBJDIR)/tmp.f90
	@if [ $(FC) = "mpiifort" ] || [ $(FC) = "mpif90" ]; then \
 	   sed "s/\!cmpi /     /g" $(OBJDIR)/tmp.f90 > $(OBJDIR)/tmp2.f90 ; \
 	   mv $(OBJDIR)/tmp2.f90 $(OBJDIR)/tmp.f90 ; \
	 fi 
	$(FC) $(OPT) -o $@ -c $(OBJDIR)/tmp.f90
	@rm -f $(OBJDIR)/tmp.f90

$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)/
	$(FC) $(OPT) -o $@ -c $<

$(OBJDIR)/:
	mkdir -p $(OBJDIR)

$(MODDIR)/:
	mkdir -p $(MODDIR)

$(EXEDIR)/:
	mkdir -p $(EXEDIR)

# Dependencies
#=============

REQMPI=
ifeq ($(FC),$(filter $(FC), mpiifort mpif90))
  REQMPI=$(OBJDIR)/module_parallelization.o
endif

$(OBJDIR)/module_mathmethods.o: $(OBJDIR)/module_constants.o

$(OBJDIR)/module_parallelization.o: $(OBJDIR)/module_mathmethods.o

$(OBJDIR)/module_nucleus.o: $(OBJDIR)/module_constants.o ${REQMPI}

$(OBJDIR)/module_basis.o: $(OBJDIR)/module_mathmethods.o $(OBJDIR)/module_nucleus.o

$(OBJDIR)/module_hamiltonian.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_wavefunctions.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_fields.o: $(OBJDIR)/module_hamiltonian.o $(OBJDIR)/module_wavefunctions.o

$(OBJDIR)/module_particlenumber.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_pairs.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_angularmomentum.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_multipoles.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_radius.o: $(OBJDIR)/module_basis.o

$(OBJDIR)/module_operators.o: $(OBJDIR)/module_particlenumber.o $(OBJDIR)/module_pairs.o \
                              $(OBJDIR)/module_angularmomentum.o $(OBJDIR)/module_multipoles.o \
                              $(OBJDIR)/module_radius.o

$(OBJDIR)/module_projection.o: $(OBJDIR)/module_fields.o $(OBJDIR)/module_operators.o 

$(OBJDIR)/module_constraints.o: $(OBJDIR)/module_nucleus.o $(OBJDIR)/module_wavefunctions.o \
                                $(OBJDIR)/module_operators.o $(OBJDIR)/module_fields.o \
                                $(OBJDIR)/module_projection.o

$(OBJDIR)/module_gradient.o: $(OBJDIR)/module_fields.o $(OBJDIR)/module_constraints.o

$(OBJDIR)/module_initialization.o: $(OBJDIR)/module_constants.o $(OBJDIR)/module_nucleus.o \
                                   $(OBJDIR)/module_hamiltonian.o $(OBJDIR)/module_wavefunctions.o \
                                   $(OBJDIR)/module_pairs.o $(OBJDIR)/module_projection.o \
                                   $(OBJDIR)/module_constraints.o $(OBJDIR)/module_gradient.o

# Debug
#======

debug:
	@echo "SRC90 = $(SRC90)"
	@echo "SRC77 = $(SRC77)"
	@echo "OBJ90 = $(OBJ90)"
	@echo "OBJ77 = $(OBJ77)"
	@echo "code  = $(code)"

# Clean up
#=========

clean:
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJDIR)/*.o  

deepclean:
	rm -rf $(MODDIR)/
	rm -rf $(OBJDIR)/

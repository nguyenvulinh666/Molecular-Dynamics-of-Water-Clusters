# Master makefile for DL_POLY_2.0 
# Author: W. Smith January 2009
# 
#=======================================================================
# Define default settings
#=======================================================================

BINROOT = ../execute
CC  = mpicc
EX = DLPOLY.X
EXE = $(BINROOT)/$(EX)
FC=undefined
SHELL=/bin/sh
TYPE=par

#=====================================================================
# Define object files

OBJ_MOD = parse_module.o setup_module.o error_module.o \
	site_module.o config_module.o pair_module.o \
        random_module.o \
	utility_module.o solvation_module.o tether_module.o \
	vdw_module.o property_module.o rigid_body_module.o \
	angles_module.o bonds_module.o shake_module.o \
	inversion_module.o dihedral_module.o core_shell_module.o \
	exclude_module.o ewald_module.o coulomb_module.o\
	external_field_module.o four_body_module.o \
	hkewald_module.o metal_module.o ensemble_tools_module.o \
	temp_scalers_module.o three_body_module.o spme_module.o \
	tersoff_module.o neu_coul_module.o \
	nlist_builders_module.o forces_module.o \
	lf_motion_module.o lf_rotation1_module.o \
	lf_rotation2_module.o vv_motion_module.o \
	vv_rotation1_module.o vv_rotation2_module.o \
	pmf_module.o integrator_module.o optimiser_module.o \
	hyper_dynamics_module.o driver_module.o \
	define_system_module.o stochastic.o \

OBJ_SRC = dlpoly.o

OBJ_PAR = basic_comms.o merge_tools.o pass_tools.o

HEAD_METAD = Plumed/metadyn.h

include Plumed/plumed.inc

OBJ_METAD = $(PLUMED_OBJECTS) images.o

#=====================================================================
# Define targets
all:
	@echo "Error - please specify a target machine!"
	@echo "Permissible targets for this Makefile are:"
	@echo "                                          "
	@echo "gfortran                   (parallel)"
	@echo "hpcx                       (parallel)"
	@echo "crayxd1                    (parallel)"
	@echo "macosx-xlf-g5-mpi          (parallel)"
	@echo "hitachi-sr2201             (parallel)"
	@echo "sg8k-mpi                   (parallel)"
	@echo "ecp-custom                 (parallel)"
	@echo "ecp-custom-opt             (parallel)"
	@echo "woodcrest                  (parallel)"
	@echo "         "
	@echo "Please examine Makefile for details"

# system specific targets follow :

#================== GNU Fortran, MPI version ==============================
gfortran:
	$(MAKE) FC="mpif90" LD="mpif90 -o" \
	LDFLAGS="-O2 -ffast-math" \
	FFLAGS="-c -O2 -ffast-math" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#================== HPCx SP Power 5 =======================================
hpcx:
	$(MAKE) FC="mpxlf" LD="mpxlf -o" \
	LDFLAGS="-O3 -q64 -qmaxmem=-1" \
	FFLAGS="-c -O3 -q64 -qmaxmem=-1 -qarch=pwr5 -qtune=pwr5 -qnosave" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#================= Cray XD1 (Portland Group) =========================
brutus:
	$(MAKE) LD="mpif90 -o" LDFLAGS="-lgoto2" \
	FC=mpif90 FFLAGS="-c -O3 " \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#============= MacOSX-XLF-G5-MPI =====================================
macosx-xlf-g5-mpi:
	$(MAKE) LD="/opt/ibmcmp/xlf/8.1/bin/xlf -o" \
	LDFLAGS="-L/opt/mpich-mx/lib -lmpich -lpmpich - L/opt/mx/lib -lmyriexpress -L/usr/lib -lSystemStubs" \
	FC="/opt/ibmcmp/xlf/8.1/bin/xlf" \
	FFLAGS="-c -qstrict -O3 -qarch=auto -qmaxmem=32768"\
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#============ Hitachi SR2201 ========================================
hitachi-sr2201:
	$(MAKE) FC=xf90 \
	FFLAGS="-c -W0,'form(fixed),opt(o(3)),langlvl(save(0))' -s,TRACE" \
	intlist.o
	$(MAKE) LDFLAGS="" LDLIBS="-lfmpi -lmpi" LD="xf90 -o" FC=xf90 \
	FFLAGS="-c -W0, 'form(fixed),opt(o(3)),langlvl(save(0))' -s,TRACE" \
	CC=xcc EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#=========== Silicon Graphics 8000 ==================================
sg8k-mpi:
	$(MAKE) LD="f90 -O3 -64 -o" FC=f90 LDFLAGS="-lmpi" \
	FFLAGS="-c -O3  -64" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#==========ECP-custom Pau for mvapich IB OPT================
ecp-custom-opt:
	$(MAKE) FC="/opt/local/mpi/mvapich-0.9.9/bin/mpif90" \
	LD="/opt/local/mpi/mvapich-0.9.9/bin/mpif90 -o" \
	LDFLAGS="-O -Mdalign" \
	FFLAGS="-c -O -Mdalign" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#===========ECP-custom Pau for mvapich IB=============================
ecp-custom:
	$(MAKE) FC="/opt/local/mpi/mvapich-0.9.9/bin/mpif90" \
	LD="/opt/local/mpi/mvapich-0.9.9/bin/mpif90 -o" \
	LDFLAGS="-O0" \
	FFLAGS="-c -O0" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#================= Woodcrest =========================================
woodcrest:
	$(MAKE) LD="mpif90 -o" LDFLAGS="" \
	FC=mpif90 FFLAGS="-c -O3" \
	EX=$(EX) BINROOT=$(BINROOT) $(TYPE)

#=====================================================================
# Default code for parallel (MPI) execution

par: check $(OBJ_MOD) $(OBJ_PAR) $(OBJ_METAD) $(OBJ_SRC)
	$(LD) $(EX) $(LDFLAGS) $(OBJ_MOD) $(OBJ_PAR) $(OBJ_METAD) $(OBJ_SRC)
	mv $(EX) $(EXE)

#=====================================================================
# Check that a machine has been specified
check:
	@if test $(FC) = "undefined";\
	then echo "You must specify a target machine!"; \
	exit 99;\
	fi

#=====================================================================
# Clean up the source directory
clean:
	rm -f $(OBJ_MOD) $(OBJ_PAR) $(OBJ_METAD) $(OBJ_SRC) *.mod

#=====================================================================
# Declare dependencies
.f.o: 
	$(FC) $(FFLAGS) $*.f

.c.o: $(HEAD_METAD)
	$(CC) -DDL_POLY -DMPI -c $(CFLAGS) $< -o $@

#=====================================================================
# Declare dependency on module files

$(OBJ_SRC): $(OBJ_MOD)

#=====================================================================
# Declare dependency for C header files

$(OBJ_METAD): $(HEAD_METAD)

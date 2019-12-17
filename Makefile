MODDIR = mod
SRCDIR = src
FC = gfortran #/usr/local/gfortran/bin/gfortran 
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS = -J $(MODDIR)
# libraries needed for linking, unused in the examples
LDFLAGS =


#DEBUG =  -mcmodel=large -march=native -ffree-line-length-none -O1 -pedantic -ffpe-trap=zero,underflow,overflow -Wall -fcheck=bounds -g
#OPTIM = -mcmodel=large -march=native -ffree-line-length-none -O3 -ffast-math -pedantic -fopenmp
DEBUG = -march=native -ffree-line-length-none -O1 -pedantic -ffpe-trap=zero,underflow,overflow -Wall -fcheck=bounds -g
OPTIM = -ffree-line-length-none -O3 -ffast-math -pedantic -fopenmp

#DEBUG = -g -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall -ffpe-trap=zero,underflow,overflow
#OPTIM = -O3 -fopenmp -funroll-loops

#-g turns on debugging
#-p turns on profiling
#-ffree-line-length-none turns off the fortran standard limit of 132 characters per line


EXE = laboetie

OBJS = \
	 $(SRCDIR)/module_precision_kinds.f90 \
	 $(SRCDIR)/module_constants.f90 \
	 $(SRCDIR)/module_input.f90 \
	 $(SRCDIR)/module_lbmodel.f90 \
	 $(SRCDIR)/module_system.f90 \
	 $(SRCDIR)/module_myallocations.f90 \
	 $(SRCDIR)/module_mathematica.f90 \
	 $(SRCDIR)/module_time.f90 \
	 $(SRCDIR)/module_io.f90 \
	 $(SRCDIR)/module_advect.f90 \
	 $(SRCDIR)/module_bounceback.f90 \
	 $(SRCDIR)/module_collision.f90 \
	 $(SRCDIR)/module_external_forces.f90 \
	 $(SRCDIR)/module_propagation.f90 \
	 $(SRCDIR)/module_transient_regime.f90 \
	 $(SRCDIR)/module_geometry.f90 \
	 $(SRCDIR)/module_moment_propagation.f90 \
	 $(SRCDIR)/module_tracers.f90 \
	 $(SRCDIR)/charges_init.f90 \
	 $(SRCDIR)/check_charge_conservation.f90 \
	 $(SRCDIR)/check_charge_distribution_equilibrium.f90 \
	 $(SRCDIR)/electrostatic_pot.f90 \
	 $(SRCDIR)/init_simu.f90 \
	 $(SRCDIR)/just_eq_smolu.f90 \
	 $(SRCDIR)/main.f90 \
	 $(SRCDIR)/poisson_nernst_planck.f90 \
	 $(SRCDIR)/smolu.f90 \
	 $(SRCDIR)/sor.f90 \
	 $(SRCDIR)/backup_phi_c_plus_c_minus.f90 \
	 $(SRCDIR)/supercell_definition.f90

# symbol '@' in front of a line makes it silent. Otherwise it is printed in terminal when called

 all: $(OBJS)
	 @mkdir -p obj mod
	 $(FC) $(FCFLAGS) $(OPTIM) -o $(EXE) $(OBJS) $(LDFLAGS)

 debug: $(OBJS)
	 @mkdir -p obj mod
	 $(FC) $(FCFLAGS) $(DEBUG) -o $(EXE) $(OBJS) $(LDFLAGS)

 clean:
	rm -vf gmon.out $(EXE) $(MODDIR)/*


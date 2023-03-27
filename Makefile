.SUFFIXES: .f90
include ./make.sys

#LFLAGS := -L$(HOME)/lib ./ARPACK_lib/libarpack.a ./ARPACK_lib/libparpack.a
#LFLAGS:=-L/usr/lib ~/scalapack-2.0.2/libscalapack.a ~/openBLAS-v0.2.13/libopenblas_nehalem-r0.2.13.a

CPP = cpp -P -C -traditional
CPPFLAGS := -DMPI

default: all

modm = mod_mpi.o mod_spring.o mod_stack.o 
modn = global_constants.o global_control.o global_context.o
plo  = plo_interface.o plo_hybf.o plo_gk.o plo_fermi.o plo_dcc.o
dmft = ctqmc_dmft.o
core = ctqmc_solver.o
lev1 = ctqmc_flavor.o ctqmc_update.o
lev2 = ctqmc_record.o ctqmc_status.o global_stream.o
lev3 = ctqmc_fourier.o ctqmc_spline.o ctqmc_util.o
dump = global_dump.o ctqmc_print.o
main = global_main.o

objects = $(modm) $(modn) $(plo) $(dmft) $(core) $(lev1) $(lev2) $(lev3) $(dump) $(main)

all: iqist_dmft

iqist_dmft: $(objects)
	  $(LINKER) $(objects) -o iqist_dmft $(LFLAGS) $(LIBS) $(LIBSMPI) $(LFLAGS)

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< $*.f90

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f iqist_dmft
	rm -f *.i90

clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat

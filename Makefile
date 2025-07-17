# This makefile is used to compile ASTR code.
# The compiler: gfortran compiler
#
#FCFLAGS=  -Wuse-without-only -g
#FC=mpif90
FC=h5pfc
#FC=ftn

SRCDIR = src user_define_module
OBJDIR = obj
BINDIR = bin
CTRDIR = cma

ifeq ($(USER),lthpc)
    export FFTW_INC_DIR=/usr/include
	export FFTW_LIB_DIR=/usr/lib/x86_64-linux-gnu
endif



FCFLAGS= -O3 -fbounds-check -I$(FFTW_INC_DIR) -lfftw3

# OPTIONS1 = -fcheck=all
OPTIONS2 = -J $(OBJDIR)
OPTIONS3 = -DHDF5
# OPTIONS4 = -DCOMB -I$(CTRDIR)/include/cantera 
# OMP = -fopenacc


EXE=Bastr2d

LIBS= -lz -lm -L$(FFTW_LIB_DIR) -lfftw3_mpi -lfftw3# -L$(CTRDIR)/lib -lcantera_fortran -lcantera -lstdc++ -pthread 
#LIBS= -lz -lm 

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

srs= commtype.F90 constdef.F90 commvar.F90 strings.F90 stlaio.F90 utility.F90  tool.F90 parallel.F90 fftwlink.F90 hdf5io.F90 readwrite.F90 solution.F90 Bastr.F90
      
OBJS=$(srs:.F90=.o)

%.o:%.F90
	@mkdir -p $(OBJDIR) 
	$(FC) $(FCFLAGS) $(INCL) $(OPTIONS1) $(OPTIONS2) $(OPTIONS3) $(OPTIONS4) $(OMP) -c -o $(OBJDIR)/$@  $<

default: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS) $(INCL) $(OMP) 

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET) $(OBJDIR)/*.mod
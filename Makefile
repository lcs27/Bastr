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


FCFLAGS= -O3 -fbounds-check -I$(FFTW_INC_DIR) -lfftw3

# OPTIONS1 = -fcheck=all
OPTIONS2 = -J $(OBJDIR)
OPTIONS3 = -DHDF5
# OPTIONS4 = -DCOMB -I$(CTRDIR)/include/cantera 
# OMP = -fopenacc


EXE=Bastr
EXE_PP=Bastrpp

LIBS= -lz -lm -L$(FFTW_LIB_DIR) -lfftw3_mpi -lfftw3# -L$(CTRDIR)/lib -lcantera_fortran -lcantera -lstdc++ -pthread 
#LIBS= -lz -lm 

TARGET = $(BINDIR)/$(EXE)
TARGET_PP = $(BINDIR)/$(EXE_PP)

VPATH = $(SRCDIR):$(OBJDIR)

srs= random.F90 cmdefne.F90 commtype.F90 constdef.F90 commvar.F90 strings.F90 stlaio.F90 utility.F90 \
	tool.F90 parallel.F90 fftwlink.F90 hdf5io.F90 readwrite.F90 solution.F90\
	Bastr.F90 Bastrpp.F90
      
OBJS=$(filter-out Bastrpp.o,$(srs:.F90=.o))
OBJS_PP=$(filter-out Bastr.o,$(srs:.F90=.o))

%.o:%.F90
	@mkdir -p $(OBJDIR) 
	$(FC) $(FCFLAGS) $(INCL) $(OPTIONS1) $(OPTIONS2) $(OPTIONS3) $(OPTIONS4) $(OMP) -c -o $(OBJDIR)/$@  $<

default: $(TARGET) $(TARGET_PP)

run: $(TARGET)

pp: $(TARGET_PP)

$(TARGET) : $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET) $(addprefix $(OBJDIR)/,$(OBJS)) $(LIBS) $(INCL) $(OMP) 

$(TARGET_PP) : $(OBJS_PP)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET_PP) $(addprefix $(OBJDIR)/,$(OBJS_PP)) $(LIBS) $(INCL) $(OMP)

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET) $(TARGET_PP) $(OBJDIR)/*.mod
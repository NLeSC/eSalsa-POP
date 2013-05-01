#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------
F77 = gfortran
F90 = gfortran
LD = gfortran
CC = cc
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes


# Adjust these to point to where netcdf is installed

# These have been loaded as a module so no values necessary
NETCDFINC = -I/var/scratch/jason/netcdf/netcdf-3.6.3-bin/include
NETCDFLIB = -L/var/scratch/jason/netcdf/netcdf-3.6.3-bin/lib

#NETCDFINC = -I/cm/shared/apps/netcdf/gcc/64/4.1.1/include
#NETCDFLIB = -L/cm/shared/apps/netcdf/gcc/64/4.1.1/lib

# Adjust these to point to where mpi is installed
#MPIINC = -I/cm/shared/apps/openmpi/gcc/64/1.4.4/include
#MPILIB = -L/cm/shared/apps/openmpi/gcc/64/1.4.4/lib64

MPIINC = -I/var/scratch/jason/OpenMPI/openmpi-1.4.2-fixed-gnu/include
MPILIB = -L/var/scratch/jason/OpenMPI/openmpi-1.4.2-fixed-gnu/lib
MPIBISLIB = -L/var/scratch/jason/climate-modelling/mpibis

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled
DHIRES               = -D_HIRES
#PRINT                = -DJASON_PRINT 
#PRINT_HALO           = -DJASON_PRINT_HALO
#PRINT_REDIST         = -DJASON_PRINT_REDIST
#PRINT_LOOP           = -DJASON_PRINT_LOOP
#TIMER                = -DJASON_TIMER 
FIX_DATA             = -DJASON_FIX_DATA
LOG_FILE             = -DJASON_SIMPLE_LOG_FILENAME
FLOW                 = -D_USE_FLOW_CONTROL
#SEND                 = -DJASON_PRINT_SEND  
FLUSH                = -DJASON_FLUSH

Cpp_opts =   \
      $(DCOUPL) $(DHIRES) $(TIMER) $(PRINT) $(PRINT_LOOP) $(LOG_FILE) $(FLOW) $(FIX_DATA) $(SEND) $(FLUSH)

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O3
# -mcmodel=medium
else
  CFLAGS := $(CFLAGS) -g -check all -ftrapuv
endif
 
#----------------------------------------------------------------------------
#
#                           FORTRAN Flags
#
#----------------------------------------------------------------------------

FBASE = $(ABI) $(NETCDFINC) $(MPIINC) $(MPI_COMPILE_FLAGS) -I$(DepDir) -fcheck-array-temporaries -Warray-temporaries 
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -O3 -fconvert=swap 
#-fmax-stack-var-size=536870912
#-mcmodel=medium
else
  FFLAGS = $(FBASE) -g -check bounds -fconvert=swap
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) -static-libgfortran
 
LIBS = $(NETCDFLIB) -lnetcdf -lcurl
 
ifeq ($(MPI),yes)
  LIBS := $(MPIBISLIB) $(MPILIB) $(MPI_LD_FLAGS) -lmpi_f77 $(LIBS) -lmpibis
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

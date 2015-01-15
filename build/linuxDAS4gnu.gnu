#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------
#
# MPILIB = -L/cm/shared/apps/openmpi/intel/64/1.4.4/lib64/
# MPIINC = -I/cm/shared/apps/openmpi/intel/64/1.4.4/include/

F77 = mpif77
F90 = mpif90
LD = mpif90 -lcurl 
CC = cc
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed
NETCDFINC = -I/cm/shared/apps/netcdf/gcc/64/4.1.1/include
NETCDFLIB = -L/cm/shared/apps/netcdf/gcc/64/4.1.1/lib

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
LOG_FILE             = -DJASON_SIMPLE_LOG_FILENAME

Cpp_opts =   \
      $(DCOUPL) $(DHIRES) $(TIMER) $(PRINT) $(PRINT_LOOP) $(LOG_FILE)

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O2
else
  CFLAGS := $(CFLAGS) -g -check all -ftrapuv
endif

CFLAGS := $(CFLAGS)
CFLAGS := $(CFLAGS) 
 
#----------------------------------------------------------------------------
#
#                           FORTRAN Flags
#
#----------------------------------------------------------------------------

FBASE = $(ABI) $(NETCDFINC) $(MPI_COMPILE_FLAGS) -I$(DepDir) 
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -O2
else
  FFLAGS = $(FBASE) -g -check bounds
endif

#DAS4 specific
FFLAGS := $(FFLAGS) -Wall -fdefault-double-8 -fdefault-real-8 -fconvert=swap -fimplicit-none -fbounds-check
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 
 
LIBS = $(NETCDFLIB) -lnetcdf 
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

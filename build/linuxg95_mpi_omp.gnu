#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------
F77 = mpif77
F90 = mpif90
LD = mpif90
CC = cc
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/awk
ABI = 
COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed

# These have been loaded as a module so no values necessary

#with -mcmodel=medium
NETCDFINC = -I/var/scratch/jason/netcdf/netcdf-4.1.1-bin-medium/include
NETCDFLIB = -L/var/scratch/jason/netcdf/netcdf-4.1.1-bin-medium/lib

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled
DHIRES                =-D_HIRES

Cpp_opts =   \
      $(DCOUPL) $(DHIRES)

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O3 -fopenmp -mcmodel=medium
else
  CFLAGS := $(CFLAGS) -g -check all -ftrapuv -fopenmp -mcmodel=medium
endif
 
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
  FFLAGS = $(FBASE) -O3 -fconvert=swap -fopenmp -mcmodel=medium
else
  FFLAGS = $(FBASE) -g -check bounds -fconvert=swap -fopenmp -mcmodel=medium
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 
 
LIBS = $(NETCDFLIB) -lnetcdf -lcurl
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) $(MPI_LD_FLAGS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS) -lgomp
 
#----------------------------------------------------------------------------

#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------

F77 = mpif90
F90 = mpif90
LD = mpif90
CC = mpicc 
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed

# These have been loaded as a module so no values necessary
NETCDFPATH = /home/ben
NETCDFINC = -I$(NETCDFPATH)/include
NETCDFLIB = -L$(NETCDFPATH)/lib

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

#DCOUPL              = -Dcoupled

Cpp_opts =   \
      $(DCOUPL)

Cpp_opts := $(Cpp_opts) -DPOSIX -DAMUSE
 
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

FFLAGS := $(FFLAGS) -Wall -fdefault-double-8 -fdefault-real-8 -fconvert=swap -fimplicit-none -fbounds-check
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 

LDFLAGS := $(LDFLAGS) $(NETCDFINC)
 
LIBS = $(NETCDFLIB) -lnetcdf -lnetcdff
 
ifeq ($(MPI),yes)
#  LIBS := $(LIBS) $(MPI_LD_FLAGS) -L/cm/shared/apps/openmpi/intel/64/1.4.4/lib64/ -lmpi 
  LIBS := $(LIBS) -L/home/ben/amuse/prerequisites -lmpi
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

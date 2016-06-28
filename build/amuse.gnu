#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------


F77 = $(MPIFC)
F90 = $(MPIFC)
LD = $(MPIFC)
CC = $(MPIFC)
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/awk
ABI = 
COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed


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

FBASE := $(ABI) $(NETCDFINC) $(MPI_COMPILE_FLAGS) -I$(DepDir) 
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

FFLAGS := $(POPFFLAGS) $(FBASE)

ifeq ($(OPTIMIZE),yes)
  FFLAGS := $(FFLAGS) -O2 
else
  FFLAGS := $(FFLAGS) -g -check bounds
endif





#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 

LDFLAGS := $(LDFLAGS) $(NETCDFINC)
 
LIBS = $(NETCDFLIB) -lnetcdf -lnetcdff
 
ifeq ($(MPI),yes)
#assuming the compiler wrapper takes care of this
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

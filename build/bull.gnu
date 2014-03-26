#-----------------------------------------------------------------------
#
# File:  bull.gnu
#
#  Contains compiler and loader options for the Bull machine 'Cartesius' at SARA in Amsterdam
#  using the intel compiler and specifies the mpi directory for communications modules.
#
#-----------------------------------------------------------------------

F77 = mpiifort -v
F90 = mpiifort -v -convert big_endian

# use gprof
#F90 = mpiifort -p -g -convert big_endian

LD = mpiifort -v
# use gprof
#LD = mpiifort -p -g

# use ita
# LD = mpiifort -g -L${VT_LIB_DIRS} -lVT ${VT_ADD_LIBS}

CC = mpiicc -v-convert big_endian
Cp = /bin/cp
Cpp = /lib/cpp -P
AWK = /usr/bin/gawk
ABI = -q64
COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.
MPI = yes

# Adjust these to point to where netcdf is installed
# These have been loaded as a module so no values necessary

NETCDFINC =  -I/hpc/sw/netcdf-4.1.3-intel-seq/include
NETCDFLIB =  -L/hpc/sw/netcdf-4.1.3-intel-seq/lib

#  Enable trapping and traceback of floating point exceptions, yes/no.
#  Note - Requires 'setenv TRAP_FPE "ALL=ABORT,TRACE"' for traceback.

TRAP_FPE = no

#------------------------------------------------------------------
#  precompiler options
#------------------------------------------------------------------

Cpp_opts = $(DCOUPL)  # Want de -P optie genereert .i files ipv .o en dat geeft een compile error

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
#CFLAGS = $(ABI) -qarch=pwr6 -qtune=pwr6 -qcache=auto
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O2
else
  CFLAGS := $(CFLAGS) -g
endif
 
#----------------------------------------------------------------------------
#
#                           FORTRAN Flags
#
#----------------------------------------------------------------------------

#FBASE = $(ABI) -qarch=pwr6 -qtune=pwr6 -qcache=auto -qzerosize $(NETCDFINC) $(MPI_COMPILE_FLAGS) -I$(DepDir) 
FBASE = $(ABI) $(NETCDFINC) $(MPI_COMPILE_FLAGS) -I$(DepDir) 
MODSUF = mod

ifeq ($(TRAP_FPE),yes)
  FBASE := $(FBASE) 
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS = $(FBASE) -O3
else
  FFLAGS = $(FBASE) -g
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) -v
 
LIBS = $(NETCDFLIB) -lnetcdf -lnetcdff
#LIBS = $(NETCDFLIB) -lnetcdf 
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) $(MPI_LD_FLAGS)
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

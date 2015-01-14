#-----------------------------------------------------------------------
#
# File:  bull.gnu
#
#  Contains compiler and loader options for the Bull machine 'Cartesius' at SARA in Amsterdam
#  using the intel compiler and specifies the mpi directory for communications modules.
#
#-----------------------------------------------------------------------



F77 = pgfortran -v -r8 -O2 -acc -Minfo=accel -ta=tesla:cuda6.0,kepler+ -byteswapio -fastsse 
F90 = pgfortran -v -r8 -O2 -acc -Minfo=accel -ta=tesla:cuda6.0,kepler+ -byteswapio -fastsse
LD = pgfortran -v -r8 -O2 -lcurl -acc -Minfo=accel -ta=tesla:cuda6.0,kepler+ -byteswapio -fastsse 
CC = pgcc -v -r8 -O2 -acc -Minfo=accel -ta=tesla:cuda6.0,kepler+ -byteswapio -fastsse
ABI = 

Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk

COMMDIR = mpi
 
#  Enable MPI library for parallel code, yes/no.
MPI = yes

# Adjust these to point to where netcdf is installed
# These have been loaded as a module so no values necessary

# FIXME?
NETCDFINC = -I/hpc/sw/netcdf-4.1.3-intel-impi-par/include
NETCDFLIB = -L/hpc/sw/netcdf-4.1.3-intel-impi-par/lib

#NETCDFINC =  -I/hpc/sw/netcdf-4.1.3-intel-seq/include
#NETCDFLIB =  -L/hpc/sw/netcdf-4.1.3-intel-seq/lib

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

CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O2
else
  CFLAGS := $(CFLAGS) -g
#  CFLAGS := $(CFLAGS) -g -check all -ftrapuv
endif

# CFLAGS := $(CFLAGS) -mcmodel=medium

 
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

# FFLAGS := $(FFLAGS) -fconvert=swap
# FFLAGS := $(FFLAGS) -mcmodel=medium
FFLAGS := $(FFLAGS) 

#----------------------------------------------------------------------------
#
#                           CUDA Flags
#
#----------------------------------------------------------------------------

# CUFLAGS = -gencode arch=compute_35,code=sm_35 -Xptxas=-v -maxrregcount=64

# CUFLAGS = -gencode arch=compute_20,code=sm_20 -Xptxas=-v

#-prec-sqrt=true -fmad=false

# ifeq ($(OPTIMIZE),yes)
#  CUFLAGS := $(CUFLAGS) 
# endif

#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) -v
 
LIBS = $(NETCDFLIB) -lnetcdff -lnetcdf -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) $(MPI_LD_FLAGS)
#  LIBS := $(LIBS) $(MPI_LD_FLAGS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#---------------------------------------------------------------------------- 

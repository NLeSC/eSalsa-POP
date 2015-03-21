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

CUDALIB = -L/cm/shared/apps/cuda55/toolkit/current/lib64/


F77 = mpif77
F90 = mpif90
LD = mpif90 -lcurl $(CUDALIB) -lcudart -lstdc++
CC = cc
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = mpi

NVCC = nvcc
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed
NETCDFINC = -I/cm/shared/apps/netcdf/gcc/64/4.1.1/include
NETCDFLIB = -L/cm/shared/apps/netcdf/gcc/64/4.1.1/lib


#with -mcmodel=medium
#NETCDFINC = -I/var/scratch/jason/netcdf/netcdf-4.1.1-bin-medium/include
#NETCDFLIB = -L/var/scratch/jason/netcdf/netcdf-4.1.1-bin-medium/lib


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

#CFLAGS := $(CFLAGS) -mcmodel=medium
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
#FFLAGS := $(FFLAGS) -mcmodel=medium

#----------------------------------------------------------------------------
#
# CUDA Flags
#
#----------------------------------------------------------------------------

#CUFLAGS = -Xptxas=-v -arch=compute_20 -code=sm_20

CUFLAGS = -gencode arch=compute_35,code=sm_35 -Xptxas=-v 

#CUFLAGS = -gencode arch=compute_20,code=sm_20 -Xptxas=-v

ifeq ($(OPTIMIZE),yes)
  CUFLAGS := -O3 $(CUFLAGS) 
endif


#  CUFLAGS := $(CUFLAGS) -Xcompiler=-mcmodel=medium


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


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
AWK = /usr/bin/gawk
ABI = 
COMMDIR = mpi
NVCC = nvcc 

#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed

# These have been loaded as a module so no values necessary
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
#FIX_DATA             = -DJASON_FIX_DATA
#LOG_FILE             = -DJASON_SIMPLE_LOG_FILENAME
FLOW                 = -D_USE_FLOW_CONTROL
#SEND                 = -DJASON_PRINT_SEND  
FLUSH                = -DJASON_FLUSH
GPU                  = -DBEN_GPU

Cpp_opts =   \
      $(DCOUPL) $(DHIRES) $(TIMER) $(PRINT) $(PRINT_LOOP) $(LOG_FILE) $(FLOW) $(FIX_DATA) $(SEND) $(FLUSH) $(PRINT_REDIST) $(GPU)

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
  CFLAGS := $(CFLAGS) -O3 -march=corei7
# -mcmodel=medium
else
  CFLAGS := $(CFLAGS) -g -check all -ftrapuv
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
  FFLAGS = $(FBASE) -O3 -march=corei7 -fconvert=swap 
#-fmax-stack-var-size=536870912
#-mcmodel=medium
else
  FFLAGS = $(FBASE) -g -check bounds -fconvert=swap
endif

#----------------------------------------------------------------------------
#
#                           CUDA Flags
#
#----------------------------------------------------------------------------

CUFLAGS = -gencode arch=compute_35,code=sm_35 -Xptxas=-v -maxrregcount=64 -gencode arch=compute_20,code=sm_20 
#CUFLAGS = -gencode arch=compute_20,code=sm_20 -Xptxas=-v

#-prec-sqrt=true -fmad=false

ifeq ($(OPTIMIZE),yes)
  CUFLAGS := $(CUFLAGS)
endif

CUFLAGS := $(CUFLAGS)
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 
 
LIBS = $(NETCDFLIB) -L/cm/shared/apps/cuda50/toolkit/current/lib64/ -lnetcdf -lcurl -lcudart -lstdc++ 
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) $(MPI_LD_FLAGS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

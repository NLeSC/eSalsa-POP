#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------

#CUDALIB = -L/cm/shared/apps/cuda40/toolkit/4.0.17/lib64/
CUDALIB = -L/cm/shared/apps/cuda50/toolkit/current/lib64/

F77 = mpif90 -O0 -g
F90 = mpif90 -O0 -g
LD = mpif90 -O0 -g  -lcurl $(CUDALIB) -lcudart -lstdc++
CC = gcc -O0 -g
Cp = /bin/cp
Cpp = cpp -P
AWK = /usr/bin/gawk
ABI = 
COMMDIR = mpi

NVCC = nvcc -O0 -g
 
#  Enable MPI library for parallel code, yes/no.

MPI = yes

# Adjust these to point to where netcdf is installed

# These have been loaded as a module so no values necessary
#NETCDFINC = -I/netcdf_include_path
#NETCDFLIB = -L/netcdf_library_path
#NETCDFINC = -I/usr/projects/climate/maltrud/local/include_coyote
#NETCDFLIB = -L/usr/projects/climate/maltrud/local/lib_coyote

#default DAS4
#NETCDFINC = -I/cm/shared/apps/netcdf/gcc/64/4.1.1/include
#NETCDFLIB = -L/cm/shared/apps/netcdf/gcc/64/4.1.1/lib

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

Cpp_opts =   \
      $(DCOUPL)

Cpp_opts := $(Cpp_opts) -DPOSIX 
 
#----------------------------------------------------------------------------
#
#                           C Flags
#
#----------------------------------------------------------------------------
 
CFLAGS = $(ABI) 

ifeq ($(OPTIMIZE),yes)
#  CFLAGS := $(CFLAGS) -O 
  CFLAGS := $(CFLAGS)  
else
  CFLAGS := $(CFLAGS) -check all -ftrapuv
endif

CFLAGS := $(CFLAGS) -mcmodel=medium
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
#  FFLAGS = $(FBASE) 
  FFLAGS = $(FBASE)
else
  FFLAGS = $(FBASE) -check bounds
endif

#DAS4 specific
FFLAGS := $(FFLAGS) -fconvert=swap
FFLAGS := $(FFLAGS) -mcmodel=medium
FFLAGS := $(FFLAGS) 
 
 
#----------------------------------------------------------------------------
#
#                           CUDA Flags
#
#----------------------------------------------------------------------------

#CUFLAGS = -Xptxas=-v -arch=compute_20 -code=sm_20

#CUFLAGS = -gencode arch=compute_35,code=sm_35 -Xptxas=-v -maxrregcount=64

CUFLAGS = -gencode arch=compute_20,code=sm_20 -Xptxas=-v

#-prec-sqrt=true -fmad=false

ifeq ($(OPTIMIZE),yes)
  CUFLAGS := $(CUFLAGS) 
endif

#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 
 
LIBS = $(NETCDFLIB) -lnetcdf
 
ifeq ($(MPI),yes)
  LIBS := $(LIBS) $(MPI_LD_FLAGS) -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

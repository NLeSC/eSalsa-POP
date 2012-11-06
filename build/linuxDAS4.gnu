#-----------------------------------------------------------------------
#
# File:  sgialtix_mpi.gnu
#
#  Contains compiler and loader options for the SGI Altix using the 
#  intel compiler and specifies the mpi directory for communications 
#  modules.
#
#-----------------------------------------------------------------------
MPILIB = -L/cm/shared/apps/openmpi/intel/64/1.4.4/lib64/
CUDALIB = -L/cm/shared/apps/cuda40/toolkit/4.0.17/lib64/

F77 = /cm/shared/apps/openmpi/intel/64/1.4.4/bin/mpif90
F90 = /cm/shared/apps/openmpi/intel/64/1.4.4/bin/mpif90
LD = /cm/shared/apps/openmpi/intel/64/1.4.4/bin/mpif90 -lcurl $(CUDALIB) -lcudart -lstdc++   -shared-intel -i-dynamic
CC = /cm/shared/apps/openmpi/intel/64/1.4.4/bin/mpicc
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
#NETCDFINC = -I/netcdf_include_path
#NETCDFLIB = -L/netcdf_library_path
#NETCDFINC = -I/usr/projects/climate/maltrud/local/include_coyote
#NETCDFLIB = -L/usr/projects/climate/maltrud/local/lib_coyote

#default DAS4
#NETCDFINC = -I/cm/shared/apps/netcdf/gcc/64/4.1.1/include
#NETCDFLIB = -L/cm/shared/apps/netcdf/gcc/64/4.1.1/lib

#intel compiler
#NETCDFINC = -I/cm/shared/apps/netcdf/intel/64/4.1.1/include
#NETCDFLIB = -L/cm/shared/apps/netcdf/intel/64/4.1.1/lib


#with -mcmodel=medium
NETCDFINC = -I/var/scratch/jason/netcdf/netcdf-4.1.1-icc-bin-medium/include
NETCDFLIB = -L/var/scratch/jason/netcdf/netcdf-4.1.1-icc-bin-medium/lib



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
  CFLAGS := $(CFLAGS) -O3
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
  FFLAGS = $(FBASE) -O3
#  FFLAGS = $(FBASE) -g
else
  FFLAGS = $(FBASE) -g -check bounds
endif

#DAS4 specific
FFLAGS := $(FFLAGS) -convert  big_endian
FFLAGS := $(FFLAGS) -mcmodel=medium -shared-intel -i-dynamic
#-i-dynamic
#FFLAGS := $(FFLAGS) 
 

#----------------------------------------------------------------------------
#
#                           CUDA Flags
#
#----------------------------------------------------------------------------

CUFLAGS = -Xptxas=-v -arch=sm_20

ifeq ($(OPTIMIZE),yes)
  CUFLAGS = $(CUFLAGS) -O3
endif
 
#----------------------------------------------------------------------------
#
#                           Loader Flags and Libraries
#
#----------------------------------------------------------------------------
 
LDFLAGS = $(ABI) 
 
LIBS = $(NETCDFLIB) -lnetcdf
 
ifeq ($(MPI),yes)
#  LIBS := $(LIBS) $(MPI_LD_FLAGS) -L/cm/shared/apps/openmpi/intel/64/1.4.4/lib64/ -lmpi 
  LIBS := $(LIBS) -L/cm/shared/apps/openmpi/intel/64/1.4.4/lib64/ -lmpi 
endif

ifeq ($(TRAP_FPE),yes)
  LIBS := $(LIBS) 
endif
 
LDLIBS = $(LIBS)
 
#----------------------------------------------------------------------------

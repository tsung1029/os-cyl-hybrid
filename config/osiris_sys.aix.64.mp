################################################################################
# System specific configuration for osiris make.
#   System:   aix 64bit
#   Compiler: xlf
#   mpi:      mp
#
# $URL:$
# $Id:$
################################################################################

BASEDIR = /ist
MPI = mp
BITS = 64


# Uncomment the following to enable parallel I/O
# PARALLEL_IO = __PARALLEL_IO__ 

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __AIX_TIMER__

# Set the following to a specific restort io backend (__RST_BINARY__, __RST_SION__)
#RST_IO = __RST_SION__
RST_IO = __RST_BINARY__



#compilers
##########
F90 = mpxlf90
cc  = xlc

MAKE = gmake
ECHO = /opt/freeware/bin/echo
SVNVERSION = /ist/subversion/current/bin/svnversion

# Fortran flags
###############

F90FLAGS = -qsuffix=cpp=f90 -q64 -qnoescape


# Production
F90FLAGS_production = $(F90FLAGS) -O3 -qarch=ppc970 -qtune=ppc970 -qcache=auto 

# Profile with xprofiler
F90FLAGS_profile    = $(F90FLAGS_production) -g -pg

# Debug
F90FLAGS_debug      = $(F90FLAGS) -g -C -qfullpath


# C flags
#########
CFLAGS = -O3 -q64


# Linker flags
##############
LDFLAGS = -q64

# Preprocessor flags
####################
FPPPREFIX = -WF,
FPPINFIX  = ,

# Link convention
#################
UNDERSCORE = FORTRANNOUNDERSCORE


# Libraries
###########
LAM_BASE =
MPI_FCOMPILEFLAGS =
MPI_FLINKFLAGS =

H5_BASE = $(BASEDIR)/hdf5/current_mp_64
H5_FCOMPILEFLAGS = -I$(H5_BASE)/lib 
H5_FLINKFLAGS = -L$(H5_BASE)/lib -L$(BASEDIR)/szip/current_64/lib -lhdf5_fortran -lhdf5 \
                -lsz -lz -lm -lgpfs

#SION_FCOMPILEFLAGS = -I
#SION_FLINKFLAGS = -L/users/u/cfpcfn/michael/packages/sion/sionlib-1.1p6/install/sionlib_aix/lib -lsion_64



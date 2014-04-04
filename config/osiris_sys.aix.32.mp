################################################################################
# System specific configuration for osiris make.
#   System:   aix 32bit
#   Compiler: xlf
#   mpi:      mp
#
# $URL:$
# $Id:$
################################################################################

BASEDIR = /ist
MPI = mp
BITS = 32

# Uncomment the following to enable parallel I/O
# PARALLEL_IO = __PARALLEL_IO__ 

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __AIX_TIMER__


#compilers
##########
F90 = mpxlf90
cc  = xlc

MAKE = gmake
ECHO = /opt/freeware/bin/echo
SVNVERSION = /ist/subversion/current/bin/svnversion

# Fortran flags
###############

F90FLAGS = -qsuffix=cpp=f90 -qnoescape


# Production
F90FLAGS_production = $(F90FLAGS) -O4 -qarch=ppc970 -qtune=ppc970 -qcache=auto 

# Profile with xprofiler
F90FLAGS_profile    = $(F90FLAGS_production) -g -pg

# Debug
F90FLAGS_debug      = $(F90FLAGS) -g -C -qfullpath -qinit=f90ptr \
                      -qflttrap=invalid:underflow:nanq:zerodivide:enable -qsigtrap \
                      -qsuppress=1518-061:1518-051 \
                      -qinitauto=FF


# C flags
#########
CFLAGS = -O3


# Linker flags
##############
#LDFLAGS = -bmaxdata:0xD0000000/dsa
LDFLAGS = -bmaxdata:0x80000000

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

H5_BASE = $(BASEDIR)/hdf5/current_mp_32
H5_FCOMPILEFLAGS = -I$(H5_BASE)/lib 
H5_FLINKFLAGS = -L$(H5_BASE)/lib -L$(BASEDIR)/szip/current_32/lib -lhdf5_fortran -lhdf5 \
                -lsz -lz -lm -lgpfs

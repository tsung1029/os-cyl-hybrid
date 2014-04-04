################################################################################
# System specific configuration for osiris make.
#   System:   MAC OS X 10.5
#   Compiler: gfortran 4.3.0 
#
# $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/config/osiris_sys.macosx.g95 $
# $Id: osiris_sys.macosx.g95 97 2007-01-03 01:34:01Z zamb $
################################################################################

MPI = ompi
BITS = 64

# Uncomment the following to enable parallel I/O
#PARALLEL_IO = __PARALLEL_IO__ 

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __MACH_TIMER__

# SIMD
# Uncomment the following line to use SIMD optimized code
SIMD = SSE

# Currently the Mac OS X assembler does not support AVX instructions. Since gcc uses the
# Mac OS X native assembler there is no simple workaround until Apple fixes this.
#SIMD = AVX

# PRECISION
#PRECISION = DOUBLE
PRECISION = SINGLE

#compilers
##########
F90 = gfortran-4.7.2
FPP = gcc -C -E -x assembler-with-cpp

cc  = gcc-4.7.2
CC  = gcc-4.7.2
#cc = icc
#CC = icc


MAKE = make
ECHO = /bin/echo
SVNVERSION = svnversion


# Fortran flags
###############

# -fno-range-check is required because of the random module. 
# gfortran has a bug that considers -2147483648 to be outside the valid
# int32 range

# -pipe makes gfortran use pipes for internal process communication (instead of files)
#       which speeds up the ocmpilation process significantly

# -ffree-line-length-none removes all constraints on line size

F90FLAGS_all = -pipe -ffree-line-length-none -fno-range-check

# OpenMP Support
F90FLAGS_all += --openmp
FPP += -D_OPENMP

# Production
# Intel Core i7 flags

F90FLAGS_production = $(F90FLAGS_all) -Ofast -march=corei7 

# Debug

# -std=f95 is too picky

F90FLAGS_debug      = $(F90FLAGS_all) -g -fbacktrace -fbounds-check \
                      -Wall -fimplicit-none -pedantic \
                      -Wimplicit-interface -Wconversion  -Wsurprising \
                      -Wunderflow -ffpe-trap=invalid,zero,overflow
                      
#-ffpe-trap=underflow this breaks vmpl

# Profile with Shark
F90FLAGS_profile    = -g $(F90FLAGS_production) 


# C flags
#########

CFLAGS_production = -Ofast -march=corei7

CFLAGS_debug = -g -Wall -pedantic -march=corei7


# Linker flags
##############

#LDFLAGS = 
LDFLAGS = -L/opt/intel/lib/intel64 -lirc


# Preprocessor flags
####################
# leave this empty to keep default.

# Link convention
#################
UNDERSCORE = FORTRANSINGLEUNDERSCORE


# Libraries
###########

MPI_FCOMPILEFLAGS = $(shell /opt/openmpi/1.6.3-gfortran/bin/mpif77 --showme:compile)
MPI_FLINKFLAGS    = $(shell /opt/openmpi/1.6.3-gfortran/bin/mpif77 --showme:link)

H5_ROOT = /opt/hdf5/1.8.10-gfortran
SZIP_ROOT = /opt/szip/2.1

H5_FCOMPILEFLAGS = -I$(H5_ROOT)/include 
H5_FLINKFLAGS    = -L$(H5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm \
                   -L$(SZIP_ROOT)/lib  -lsz

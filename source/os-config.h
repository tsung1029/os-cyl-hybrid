! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define __HAS_MPI_IN_PLACE__

! Use OpenMP 
!#define _OPENMP 

! Use parallel I/O
!#define __PARALLEL_IO__ 

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
!#define SIMD_AVX
!#define SIMD_BGP
!#define SIMD_BGQ

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not support use of namelist I/O with internal files (e.g. Pathscale 3.2)
! which is a Fortran 2003 feature
!#define __NO_NML_INTFILE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------

#ifndef __DISTRO__

! Include Ionization module
#define __HAS_IONIZATION__

! Include binary collisions module
#define __HAS_COLLISIONS__

! Include particle tracking module
#define __HAS_TRACKS__

! Include perfectly matched layers boundary conditions for EMF
#define __HAS_PML__

#endif


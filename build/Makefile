################################################################################
# Multiplatform Makefile for osiris 2.0 source files
#
# $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/config/Makefile.build $
# $Id: Makefile.build 555 2013-04-02 13:06:49Z zamb $
################################################################################

         
# Object Files
##############
OBJSBASE = os-sys-multi-c.o os-sys-multi.o memory.o os-param.o \
           os-units.o os-math.o os-files.o 

OBJSLOG = os-logprof.o

OBJSALL =         os-stringutil.o os-util.o  os-restart.o\
  os-random.o os-nconf.o    os-fparser.o  \
  os-grid-define.o os-grid-memory.o os-grid-parallel.o os-grid-report.o os-grid.o \
  os-space.o    os-tstep.o  os-dutil-hdf5.o os-dutil.o  \
  os-vdf-define.o os-vdf-memory.o os-vdf-reportfile.o os-vdf-comm.o os-vdf-math.o os-vdf-smooth.o \
  os-vdf-average.o os-vdf.o os-vdf-report.o os-vdf-interpolate.o \
  os-wall-define.o os-wall-comm.o os-wall.o  os-time.o     os-timeavg.o   \
  os-cyl-modes.o os-current-define.o \
  os-emf-define.o os-emf-lindman.o os-emf-vpml.o os-emf-boundary.o \
  os-emf-marder.o os-emf-pgc.o os-emf-cyl-modes.o os-emf-solver.o os-emf-interpolate.o \
  os-emf-interp-cm.o os-emf-poynting.o os-emf-psi.o os-emf-diagnostics.o \
  os-emf-gridval.o os-emf.o \
  os-current-diagnostics.o  os-current.o \
  os-spec-define.o os-spec-memory.o os-spec-v3dep.o os-spec-v3int.o \
  os-spec-push-c.o os-spec-profile.o os-spec-tracks.o os-spec-udist.o os-spec-util.o \
  os-spec-comm.o  os-spec-current-c.o os-spec-current.o os-spec-charge.o \
  os-spec-piston.o  os-spec-rawdiag.o os-spec-phasespace.o os-spec-cyl-modes.o os-spec-boundary.o \
  os-spec-loadbalance.o  os-spec-pgc.o os-spec-radcool.o  os-spec-diagnostics.o  os-species.o   \
  os-cathode-define.o os-cathode.o os-cross.o os-dneutral.o \
  os-neutral-define.o  os-neutral.o \
  os-spec-collisions-define.o os-spec-collisions.o \
  os-particles-define.o os-particles-charge.o os-particles.o\
  os-zpulse.o  os-antenna.o \
  os-antenna_array.o
  
OBJSMAIN =  os-main.o

OBJS    = $(OBJSBASE) $(OBJSLOG) $(OBJSALL) $(OBJSMAIN)


# ioperf tool
OBJS_IOPERF = $(OBJSBASE) $(OBJSLOG) $(OBJSALL) ioperf.o


# shell
#######
SHELL = /bin/sh
PYTHON = python

############################
## configuration file     ##
############################
-include ../config/osiris_config 

# default preprocessor
ifeq ("$(FPP)","")
  FPP = gcc -C -E -x assembler-with-cpp
endif


############################
# build dependencies      ##
############################
-include f_deps.mk

.PHONY: depend
depend:
	@if [ "$(ENABLE_DEPS)" = "1" ]; then \
          ../tools/findFDep.sh "$(MOD_NAME)" "$(FPP)"; \
          else echo "" > f_deps.mk; fi


# get osiris revision
#####################
ifeq ("$(VERSION)","")
VERSION := r$(shell svnversion -c ../source/ | awk -F: '{ print $$2 }' )
endif


# Compiler flags
ifneq ("$(CFLAGS_$(COMPILATION_TYPE))","")
  # New style, accepts optimization flags for C code
  CF   = $(CFLAGS_$(COMPILATION_TYPE))
else
  # Old style, uses the same C flags for all targets
  CF   = $(CFLAGS)
endif

# Add Fortran underscore definition
CF += -D$(UNDERSCORE)

F90F = $(F90FLAGS_$(COMPILATION_TYPE)) $(MPI_FCOMPILEFLAGS) $(H5_FCOMPILEFLAGS)

FPPF = -DP_X_DIM=$(OS_DIM) -DOS_REV=\"$(VERSION)\" -D$(UNDERSCORE)
LDF  = $(LDFLAGS) $(H5_FLINKFLAGS) $(MPI_FLINKFLAGS) $(SION_FLINKFLAGS)


# Precision flags
ifneq ("$(PRECISION)","")
  ifneq ("$(PRECISION)","SINGLE")
   ifneq ("$(PRECISION)","DOUBLE")
     $(error When set, PRECISION must be either SINGLE or DOUBLE. Please correct your configuration file)
   endif  
  endif  
  CF   +=-DPRECISION_$(PRECISION)
  FPPF +=-DPRECISION_$(PRECISION)
endif


# SIMD flags
ifneq ("$(SIMD)","")
  ifeq ("$(PRECISION)","")
    $(error SIMD code selected, but PRECISION was not set. Please correct your configuration file)
  endif
  CF   +=-DSIMD -DSIMD_$(SIMD)
  FPPF +=-DSIMD -DSIMD_$(SIMD)
endif

# Parallel IO/Flags
ifneq ("$(PARALLEL_IO)","")
  FPPF+=-D$(PARALLEL_IO)
endif

# Allocatable array extensions
ifneq ("$(ALLOC_EXT)","")
  FPPF+=-D$(ALLOC_EXT)
endif

# System timers
ifneq ("$(TIMER)","")
  CF+=-D$(TIMER)
endif

# Disable file flushing
ifneq ("$(NO_FILE_FLUSH)","")
  FPPF+=-D$(NO_FILE_FLUSH)
endif

# Enable log files
ifneq ("$(USE_LOG)","")
  FPPF+=-D$(USE_LOG)
endif

# Restart io backend
ifneq ("$(RST_IO)","")
  FPPF+=-D$(RST_IO)
endif

# rules
########

#c
%.o : ../source/%.c ../config/osiris_config
	$(cc) $(CF) $(INCPATH) -c $<

#cpp
%.o : ../source/%.cpp ../config/osiris_config
	$(CC) $(CF) $(INCPATH) -c $<

#f90 without and with preprocessor
%.o :../source/%.f90 ../config/osiris_config ../source/os-param.f90
	$(FPP) $(FPPF) -I../source $< > $(<F)
	$(F90) $(F90F) $(INCPATH) -c $(<F)

osiris.e : $(OBJS) ../config/osiris_config
	$(F90)  $(F90F) -o $@ $(OBJS) $(LDF)
	chmod 755 osiris.e
	cp osiris.e ../bin/osiris-$(VERSION)-$(OS_DIM)D.e
	ln -f -s osiris-$(VERSION)-$(OS_DIM)D.e ../bin/osiris-$(OS_DIM)D.e

.PHONY: clean
clean :
	touch a.o && rm -f *.o
	touch a.f90 && rm -f *.f90

ioperf : $(OBJS_IOPERF)
	$(F90)  $(F90F) -o $@ $(OBJS_IOPERF) $(LDF)
	chmod 755 ioperf
	cp ioperf ../bin/ioperf-$(BITS)-$(MPI)-$(COMPILATION_TYPE)-$(OS_VERSION)-$(OS_DIM)D
	ln -f -s ioperf-$(BITS)-$(MPI)-$(COMPILATION_TYPE)-$(OS_VERSION)-$(OS_DIM)D ../bin/ioperf-$(OS_DIM)D

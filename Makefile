################################################################################
# Multiplatform Makefile for osiris 2.0 distribution (main makefile)
#
# $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/Makefile $
# $Id: Makefile 337 2010-03-26 14:27:00Z michael $
################################################################################

# globals
#########
SHELL = /bin/sh

# targets
#########

all: ./config/osiris_config
	cd build && $(MAKE) osiris.e

1d: ./config/osiris_config
	./configure -d 1 && cd build && $(MAKE) osiris.e

2d: ./config/osiris_config
	./configure -d 2 && cd build && $(MAKE) osiris.e

3d: ./config/osiris_config
	./configure -d 3 && cd build && $(MAKE) osiris.e

dist: ./config/osiris_config 1d 2d 3d

distall: ./config/osiris_config
	./configure -t debug && $(MAKE) dist
	./configure -t production && $(MAKE) dist

production: ./config/osiris_config
	./configure -t production && cd build && $(MAKE) osiris.e

profile: ./config/osiris_config
	./configure -t profile && cd build && $(MAKE) osiris.e

debug: ./config/osiris_config
	./configure -t debug && cd build && $(MAKE) osiris.e

ioperf: ./config/osiris_config
	cd build && $(MAKE) ioperf

# configuration
###############

config/osiris_config:
	./configure -l


# clean etc.
############

.PHONY: clean
clean: 
	cd build && touch osiris.e && rm -f osiris.e
	cd build && touch osiris.e.dSYM && rm -rf osiris.e.dSYM
	cd build && $(MAKE) clean
	cd build && touch mod.mod && rm -f *.mod
	cd build && touch ioperf && rm -f ioperf

.PHONY: clean-exec
clean-exec: 
	rm -f build/osiris.e bin/osiris.1D.e bin/osiris.2D.e bin/osiris.3D.e

.PHONY: clean-objs
clean-objs:
	cd build && $(MAKE) clean
	cd build && rm -f *.mod

.PHONY: distclean
distclean:
	rm -rf build
	rm -rf bin
	rm -f config/osiris_config

.PHONY: test
test: ./bin/osiris-1D.e ./bin/osiris-2D.e ./bin/osiris-3D.e
	mkdir -p test
	cd test && ln -fs ../config/Makefile.test Makefile
	cd test && $(MAKE) test
	rm -rf test

./bin/osiris-1D.e :
	./configure -d 1 && cd build && $(MAKE) osiris.e

./bin/osiris-2D.e :
	./configure -d 2 && cd build && $(MAKE) osiris.e

./bin/osiris-3D.e :
	./configure -d 3 && cd build && $(MAKE) osiris.e

.PHONY: help
help:
	@echo "NOT UP TO DATE::::::::::::::::::"
	@echo "Osiris Makefile options"
	@echo
	@echo "all          - (default) Build build/osiris.e, bin/osiris.?D.e and targets"
	@echo "1d           - Build build/osiris.e and bin/osiris.1D.e"
	@echo "2d           - Build build/osiris.e and bin/osiris.2D.e"
	@echo "3d           - Build build/osiris.e and bin/osiris.3D.e"
	@echo "dist         - Build bin/osiris.1D.e bin/osiris.2D.e and bin/osiris.3D.e"
	@echo "clean        - Remove executables and object files"
	@echo "clean-exec   - Remove executables"
	@echo "clean-objs   - Remove object files"
	@echo "help         - Display this message"

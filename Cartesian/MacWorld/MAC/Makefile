#
# Makefile for MAC platform
#
LICENSE :=
CCC := gcc
OPT_TEST := opt0
RM    := rm -fr
MKDIR := mkdir -p
DATE  := date

# Get the architecture from etc/xxx
SYS:=$(shell perl -I tools/Perl tools/Perl/arch.pl $(CCC))
ifndef SYS
$(shell perl -I tools/Perl tools/Perl/arch.pl -v $(CCC))
$(error Unable to find a valid architecture for this computer)
endif

CONFIG=bin/init
BEGIN_TARGET=@printf "****************************************\n* \
Begin target : %s \n****\n" 
END_TARGET=@printf "****\n* \
target %s completed\n****************************************\n"

# Get extension for dynamic libraries from etc/xxx
DYNAMIC_LIB_EXT := $(shell perl -I tools/Perl tools/Perl/arch.pl $(CCC) -getvariable DYNAMIC_LIB_EXT)

libdir = lib/$(SYS)
abslibdir := $(shell perl -MFile::Spec -e 'print File::Spec->rel2abs("$(libdir)")')
lib0 = $(libdir)/libmac0$(DYNAMIC_LIB_EXT)
lib1 = $(libdir)/libmac1$(DYNAMIC_LIB_EXT)
lib2 = $(libdir)/libmac2$(DYNAMIC_LIB_EXT)
libg = $(libdir)/libmacg$(DYNAMIC_LIB_EXT)


##########################
# MACRO TARGETS          #
##########################

default:help


##################################
# GNU-like MACRO TARGETS         #
##################################

all:libraries

clean:cleanlib

cpinclude:cpinc


#######################
# Libraries builder   #
#######################

lib0 : $(lib0)
lib1 : $(lib1)
lib2 : $(lib2)
libg : $(libg)
#libpg : $(libdir)/libmacpg$(DYNAMIC_LIB_EXT)

libraries : lib1 lib0
	$(END_TARGET) $@

extra-libraries : lib2 libg
	$(END_TARGET) $@

.PRECIOUS:$(libdir)/%/lib$(DYNAMIC_LIB_EXT) $(libdir)/%/Makefile

#$(libdir)/libmacg$(DYNAMIC_LIB_EXT): $(libdir)/dbg/Makefile FORCE
#	mac build -make $(MAKE) -archive ../$(@F) $(<D)

$(libdir)/libmac%$(DYNAMIC_LIB_EXT): $(libdir)/opt%/Makefile FORCE
	mac build -make $(MAKE) -archive ../$(@F) $(<D)

$(libdir)/%/Makefile: 
	$(MKDIR) $(@D)
	mac depend -verbose -mMAC -compiler $(CCC) $* $(@D)


#######################
# cleaner             #
#######################

cleanlib: 
	$(RM) lib


################
# Copy include #
################

cpinc:
	mac cpinclude -verbose -mMAC -mAPI


########
# Help #
########

help:
	@echo
	@echo "Usage: $(MAKE) CCC=<arg> <target>"
	@echo "<target>"
	@echo "    all            : build all libraries"
	@echo "    lib0           : build library lib0"
	@echo "    lib1           : build library lib1"	
	@echo "    lib2           : build library lib2"	
	@echo "    libg           : build library lib0"	
	@echo "    libraries      : build libraries lib0 and lib1"
	@echo "    extra-libraries: build libraries lib3 and libg"
	@echo "    clean          : remove all generated files and libraries"
	@echo "    cpinclude      : copy include files in include/"
	@echo "<arg> : compiler system (default: gcc)"
	@echo

FORCE:
# DO NOT DELETE

##################################################################
#        Standard GNU Makefile for MAC Applications              #
##################################################################


# Sets here the directories where the mac scripts will recursively
# find the source files of the application
SRC = RegressionTests Misc Solvers FluidSolidInterface PostProcessing

# default compilation system (used to find the architecture)
CCC := gcc

# Test if the current architecture is valid for MAC
SYS    := $(shell mac arch $(CCC))

ifndef SYS
$(shell mac arch -v $(CCC))
$(error Unable to find a valid architecture for this computer)
endif

# directory where objects files and binaries are stored
libdir := lib/$(SYS)

# targets that are not the name of a file. There are just names
# for some commands to be executed on an explicit request
.PHONY : exe all clean exe% FORCE

# Tells make not to remove intermediate target such as Makefile(s) in $(libdir)
.SECONDARY :

# the default target
all : help

clean:
	$(BEGIN_TARGET) clean
	$(RM) lib
	$(END_TARGET) clean


##################################################################
# Executables                                                    #
# the user can make predefined targets :                         #
#       make exe2 => application opt2, mac opt1                  #
#       make exe0 => application opt0, mac opt0                  #
# the user may also build non predefined targets :               #
#       make exe2 MACLIB=mac0 => application opt2, mac opt0      #
#       make exeg MACLIB=mac2 => application dbg, mac opt2       #
# so more generally,                                             #
#       make exeOOO MACLIB=PPP => application OOO, mac PPP       #
# some very special targets are :                                #
#       make exegg => application dbg, mac dbg                   #
#       make exepg => application optpg, mac mac0                #
#       make exepgpg => application optpg, mac optpg             #
#            optpg option is used to profile an application      #
##################################################################

# the application is linked with the MACLIB version of the MAC librarie
MACLIB = mac0

# default target generates 3 executables
exe: exe0 exe2 exeg

# overrides MACLIB for some known targets
exe2    : MACLIB=mac1
exeg    : MACLIB=mac1
exegg   : MACLIB=macg
exepgpg : MACLIB=$(MAC_HOME)/$(libdir)/optpg

$(libdir)/optg/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	mac depend -l mac1 -compiler $(CCC) dbg $(@D) $(SRC)

$(libdir)/optg%/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	mac depend -l mac$* -compiler $(CCC) dbg $(@D) $(SRC)

$(libdir)/optpg%/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	mac depend -l $(MACLIB) -compiler $(CCC) optpg $(@D) $(SRC)

$(libdir)/%/Makefile:
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	mac depend -verbose -l $(MACLIB) \
	-I $(XERCESCROOT)/include \
	-I $(MACWORLD_MPI_INCDIR) \
	-I $(BASILISK_DLMFD) \
	-I $(GRAINS_HOME)/Grains/inc \
	-path $(GRAINS_HOME)/Grains/lib$(GRAINS_FULL_EXT) -l grains \
	-path $(XERCESCROOT)/lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION) -l xerces-c -l xerces-depdom \
	-path $(MACWORLD_GFORTRAN_LIBDIR) -l $(MACWORLD_GFORTRAN_LIBS) \
	-compiler $(CCC) $* $(@D) $(SRC)
	$(END_TARGET) $@
	
exe% : $(libdir)/exe%
	$(END_TARGET) $@

$(libdir)/exe% : $(libdir)/opt%/Makefile FORCE
	$(BEGIN_TARGET) $@
	mac build -exe $(libdir)/opt$*
	if [ ! -f $@ ] ; then cd $(@D) ; ln -s opt$*/exe $(@F) ; fi
	$(END_TARGET) $@


##################################################################
# Help                                                           #
##################################################################

ifdef CALLER
COMMANDE := $(CALLER)
else
COMMANDE := $(MAKE)
endif

help:
	@echo "---------------------------------------------------------"
	@echo "Usage : $(COMMANDE) CCC=<arg> <target>"
	@echo "---------------------------------------------------------"
	@echo " <arg>"
	@echo "   compiler system (default gcc)"
	@echo "---------------------------------------------------------"
	@echo " <target>"
	@echo "   clean : remove all generated files"
	@echo "   exe0  : build executable : opt0  ; libmac0.so"
	@echo "   exe2  : build executable : opt2  ; libmac1.so"
	@echo "   exeg  : build executable : dbg   ; libmac1.so"
	@echo "   exegg : build executable : dbg   ; libmacg.so"
	@echo "   exepg : build executable : optpg ; libmac0.so"
	@echo "   exepgpg : build full profiled executable"
	@echo "   exe   : same as exe0, exe2 and exeg"
	@echo "---------------------------------------------------------"
	@echo " architecture : " $(SYS)
	@echo "---------------------------------------------------------"

FORCE:

##################################################################
# internal commands                                              #
##################################################################

RM    := rm -fr
MKDIR := mkdir -p
BEGIN_TARGET := @printf "****************************************\n* \
Begin target : %s \n****\n" 
END_TARGET := @printf "****\n* \
target %s completed\n****************************************\n"


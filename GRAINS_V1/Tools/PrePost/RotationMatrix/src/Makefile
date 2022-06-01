# Target directories
EXEDIR := ../bin$(GRAINS_FULL_EXT)
OBJDIR := ../obj$(GRAINS_FULL_EXT)


# Source and object files
MAINSRC := rotationmatrix.cpp
MAINOBJ := $(MAINSRC:%.cpp=$(OBJDIR)/%.o)


# Dependencies of
GRAINSINC   := $(GRAINS_ROOT)/include
GRAINSLIB   := $(GRAINS_ROOT)/lib$(GRAINS_FULL_EXT)
GRAINSFLAGS := -L$(GRAINSLIB) -lgrains

XERCESINC   := $(XERCESC_ROOT)/include
XERCESLIB   := $(XERCESC_ROOT)/lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
XERCESFLAGS := -L$(XERCESLIB) -lxerces-c -lxerces-depdom


# Compiler
COMP  := $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
LINK  := $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
DEBUG := -g


# Compiler flags
INCFLAGS  := $(XERCESINC:%=-I%) $(GRAINSINC:%=-I%)
COMPFLAGS := $(DEBUG) $(INCFLAGS) -fPIC
LINKFLAGS := $(DEBUG)


# -----------------------------------------------------------------------------
# Production rules
.SILENT:

all : MESS $(EXEDIR)/rotationmatrix
	@echo 'PrePost tool rotationmatrix built !'
	@echo

$(EXEDIR)/rotationmatrix : $(MAINOBJ)
	@echo 'Building PrePost tool rotationmatrix'
	$(LINK) $(MAINOBJ) -o $@ $(LINKFLAGS) $(GRAINSFLAGS) $(XERCESFLAGS)

MESS   :
	@echo '***********************************************'
	@echo 'Compiling rotationmatrix in' $(GRAINS_BITS_DEFAULT) 'bits with' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
	@echo '***********************************************'

$(OBJDIR)/%.o : %.cpp
	@echo 'Compiling '$<
	$(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@


# DEPENDENCIES
UTLDIR  := $(GRAINS_HOME)/Tools/Perl
INCLUDE := $(GRAINSINC) $(XERCESINC) $(GRAINS_SYSINC)
DEPEND  := $(INCLUDE:%=-I%)

depend :
	@echo 'Generating dependencies'
	touch ./Maketmp
	makedepend -f ./Maketmp $(DEPEND) $(MAINSRC)
	$(UTLDIR)/makeadd.pl $(OBJDIR)/ ./Maketmp ./Makedepend
	rm -f ./Maketmp ./Maketmp.bak

include ./Makedepend      # DO NOT DELETE

# CLEAN
clean :
	rm -f $(MAINOBJ)
	rm -f $(EXEDIR)/*
  

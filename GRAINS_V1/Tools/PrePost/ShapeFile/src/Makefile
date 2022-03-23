# Target directories
EXEDIR := ../bin$(GRAINS_FULL_EXT)
OBJDIR := ../obj$(GRAINS_FULL_EXT)


# Source and object files
MAINSRC := shapefile.cpp
MAINOBJ := $(MAINSRC:%.cpp=$(OBJDIR)/%.o)


# Compiler
COMP  := $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
LINK  := $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
DEBUG := -g


# Compiler flags
INCFLAGS  := 
COMPFLAGS := $(DEBUG) $(INCFLAGS) -fPIC
LINKFLAGS := $(DEBUG)


# -----------------------------------------------------------------------------
# Production rules
.SILENT:

all : MESS $(EXEDIR)/shapefile
	@echo 'PrePost tool shapefile built !'
	@echo

$(EXEDIR)/shapefile : $(MAINOBJ)
	@echo 'Building PrePost tool shapefile'
	$(LINK) $(LINKFLAGS) $(MAINOBJ) -o $@

MESS   :
	@echo '*************************************************'
	@echo 'Compiling shapefile in' $(GRAINS_BITS_DEFAULT) 'bits with' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
	@echo '*************************************************'

$(OBJDIR)/%.o : %.cpp
	@echo 'Compiling '$<
	$(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@


# DEPENDENCIES
UTLDIR  := $(GRAINS_HOME)/Tools/Perl
INCLUDE := $(GRAINS_SYSINC)
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
  

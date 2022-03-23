# Target directories
EXEDIR := ../bin$(GRAINS_FULL_EXT)
OBJDIR := ../obj$(GRAINS_FULL_EXT)


# Source and object files
MAINSRC := main.cpp
MAINOBJ := $(MAINSRC:%.cpp=$(OBJDIR)/%.o)


# Include and library from Grains
GRAINSINC   := $(GRAINS_ROOT)/include
GRAINSLIB   := $(GRAINS_ROOT)/lib$(GRAINS_FULL_EXT)
GRAINSFLAGS := -L$(GRAINSLIB) -lgrains


# Compiler
# Comment: the general CC flags and MPI flags are defined in the 
# environment file and should not be altered unless you are a hard core
# Grains developer
COMP  := $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) $(GRAINS_GENCCFLAGS) $(GRAINS_MPICCFLAGS)
LINK  := $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) $(GRAINS_GENCCFLAGS) $(GRAINS_MPICCFLAGS)
DEBUG  := -g


# XML 
MACRO       = -D_XML
XERCESINC   := $(XERCESC_ROOT)/include
XERCESLIB   := $(XERCESC_ROOT)/lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
XERCESFLAGS := -L$(XERCESLIB) -lxerces-c -lxerces-depdom


# MPI
MPIINC	    := $(GRAINS_MPI_INCDIR)
MPILIB      := $(GRAINS_MPI_LIBDIR)
MPIFLAGS    := -L$(MPILIB) $(LIBMPI_FOR_GRAINS___)


# Other macros
MACRO      += -D_TIMER


# Compiler flags
INCFLAGS  := $(XERCESINC:%=-I%) $(GRAINSINC:%=-I%) $(MPIINC:%=-I%)
COMPFLAGS := $(DEBUG) $(INCFLAGS)
LINKFLAGS := $(DEBUG)


# -----------------------------------------------------------------------------
# Production rules
.SILENT:

all : grains

# EXECUTABLE
grains : MESS START $(EXEDIR)/grains
	@echo '----------------------'
	@echo 'Software grains built!'
	@echo '----------------------'

$(EXEDIR)/grains : app
	@echo 'Building software grains'
	$(LINK) $(LINKFLAGS) $(XERCESFLAGS) $(GRAINSFLAGS) $(MPIFLAGS) $(MAINOBJ) -o $@

MESS   :
	@echo '***********************************************************'
	@echo 'Binary compiled in' $(GRAINS_BITS_DEFAULT) 'bits with' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION) 'and' $(GRAINS_MPI_DISTRIB)-$(GRAINS_MPI_VERSION)
	@echo '***********************************************************'

START  :
	@echo 'MPI libraries = '$(LIBMPI_FOR_GRAINS___)
	@echo 'Compilation flags = $(GRAINS_GENCCFLAGS)$(GRAINS_MPICCFLAGS)$(COMPFLAGS)'
	@echo 'Macros = $(MACRO)'
	@echo
	
app    :  Grains $(MAINOBJ)

Grains :
	( cd ../../Grains; make; cd ../Main/src )

$(OBJDIR)/%.o : %.cpp
	@echo 'Compiling '$<
	$(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@


# DEPENDENCIES
UTLDIR  := $(GRAINS_HOME)/Tools/Perl
DEPEND  := -I$(GRAINSINC) $(GRAINS_SYSINC:%=-I%)

depend :
	@echo 'Generating dependencies'
	touch ./Maketmp
	makedepend -f ./Maketmp $(DEPEND) $(MAINSRC)
	$(UTLDIR)/makeadd.pl $(OBJDIR)/ ./Maketmp ./Makedepend
	rm -f ./Maketmp ./Maketmp.bak

include ./Makedepend      # DO NOT DELETE

# CLEAN
clean : 
	( cd ../../Grains; make clean; cd ../Main/src )
	@echo 'Main   : deleting binary and object files'
	rm -f $(MAINOBJ)
	rm -f $(EXEDIR)/*	

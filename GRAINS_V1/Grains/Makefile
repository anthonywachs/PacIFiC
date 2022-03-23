# Source and include directories
INCDIR =
INCDIROUT := include
ALLINC = 

SRCDIRBASE := Base/src
INCDIRBASE := Base/include
INCDIR += $(INCDIRBASE)
ALLINC += $(wildcard $(INCDIRBASE)/*.hh)

SRCDIRCOMPO := Component/src
INCDIRCOMPO := Component/include
INCDIR += $(INCDIRCOMPO)
ALLINC += $(wildcard $(INCDIRCOMPO)/*.hh)

SRCDIRFLUID := Fluid/src
INCDIRFLUID := Fluid/include
INCDIR += $(INCDIRFLUID)
ALLINC += $(wildcard $(INCDIRFLUID)/*.hh)

SRCDIRGEOM := Geometry/src
INCDIRGEOM := Geometry/include
INCDIR += $(INCDIRGEOM)
ALLINC += $(wildcard $(INCDIRGEOM)/*.hh)

SRCDIRPP := PostProcessing/src
INCDIRPP := PostProcessing/include
INCDIR += $(INCDIRPP)
ALLINC += $(wildcard $(INCDIRPP)/*.hh)

SRCDIRTIME := TimeIntegration/src
INCDIRTIME := TimeIntegration/include
INCDIR += $(INCDIRTIME)
ALLINC += $(wildcard $(INCDIRTIME)/*.hh)

SRCDIRXML := XML/src
INCDIRXML := XML/include
INCDIR += $(INCDIRXML)
ALLINC += $(wildcard $(INCDIRXML)/*.hh)


# Target directories
OBJDIR := obj$(GRAINS_FULL_EXT)
LIBDIR := lib$(GRAINS_FULL_EXT)


# Source files per package
# _WP_ means "with path"
SRCFILES_WP_BASE := $(wildcard $(SRCDIRBASE)/*.cpp)
SRCFILESBASE := $(subst $(SRCDIRBASE)/,,$(SRCFILES_WP_BASE))

SRCFILES_WP_COMPO := $(wildcard $(SRCDIRCOMPO)/*.cpp)
SRCFILESCOMPO := $(subst $(SRCDIRCOMPO)/,,$(SRCFILES_WP_COMPO))

SRCFILES_WP_FLUID := $(wildcard $(SRCDIRFLUID)/*.cpp)
SRCFILESFLUID := $(subst $(SRCDIRFLUID)/,,$(SRCFILES_WP_FLUID))

SRCFILES_WP_GEOM := $(wildcard $(SRCDIRGEOM)/*.cpp)
SRCFILESGEOM := $(subst $(SRCDIRGEOM)/,,$(SRCFILES_WP_GEOM))

SRCFILES_WP_PP := $(wildcard $(SRCDIRPP)/*.cpp)
SRCFILESPP := $(subst $(SRCDIRPP)/,,$(SRCFILES_WP_PP))

SRCFILES_WP_TIME := $(wildcard $(SRCDIRTIME)/*.cpp)
SRCFILESTIME := $(subst $(SRCDIRTIME)/,,$(SRCFILES_WP_TIME))

SRCFILES_WP_XML := $(wildcard $(SRCDIRXML)/*.cpp)
SRCFILESXML := $(subst $(SRCDIRXML)/,,$(SRCFILES_WP_XML))


# Substitutions for objects
# Here we use the base file names combined to VPATH
OBJS_BASE := $(SRCFILESBASE:%.cpp=$(OBJDIR)/%.o)
OBJS_COMPO := $(SRCFILESCOMPO:%.cpp=$(OBJDIR)/%.o)
OBJS_FLUID := $(SRCFILESFLUID:%.cpp=$(OBJDIR)/%.o)
OBJS_GEOM := $(SRCFILESGEOM:%.cpp=$(OBJDIR)/%.o)
OBJS_PP := $(SRCFILESPP:%.cpp=$(OBJDIR)/%.o)
OBJS_TIME := $(SRCFILESTIME:%.cpp=$(OBJDIR)/%.o)
OBJS_XML := $(SRCFILESXML:%.cpp=$(OBJDIR)/%.o)


# Defining VPATH for automatic search of source files in the 
# different source directories 
# This enables to have a single rule/target for all source files in 
# various directories and all objects in a single directory
VPATH := $(SRCDIRBASE):$(SRCDIRCOMPO):$(SRCDIRFLUID):$(SRCDIRGEOM):$(SRCDIRPP):$(SRCDIRTIME):$(SRCDIRXML)


# Compiler
# Comment: the general CC flags and MPI flags are defined in the 
# environment file and should not be altered unless you are a hard core
# Grains developer
COMP  := $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) $(GRAINS_GENCCFLAGS) $(GRAINS_MPICCFLAGS)
LINK  := $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) $(GRAINS_GENCCFLAGS) $(GRAINS_MPICCFLAGS)
DEBUG  := -g


# XML 
MACRO       := -D_XML
XERCESINC   := $(XERCESC_ROOT)/include


# MPI
MPIINC	    := $(GRAINS_MPI_INCDIR)
MPILIB      := $(GRAINS_MPI_LIBDIR)
#MPIFLAGS    = -L$(MPILIB) $(LIBMPI_FOR_GRAINS___)
export LIBMPI_FOR_GRAINS___ = `echo ${GRAINS_MPI_LIBS} | sed 's%[^ ][^ ]*%-l&%g'`
MPIFLAGS    := -L$(MPILIB) $(LIBMPI_FOR_GRAINS___)


# Zlib
ZLIB := /usr/lib$(GRAINS_BITS_EXT)
ZFLAGS := -L$(ZLIB) -lz


# Compiler flags
INCFLAGS  := $(INCDIR:%=-I%) $(XERCESINC:%=-I%) $(MPIINC:%=-I%)
COMPFLAGS := $(DEBUG) $(INCFLAGS)



#------------------------------------------------------------------------------
# Production rules
.SILENT:

all : grains
	@echo '--------------------------------'

grains : MESS CPINC START XML GEOM COMPO BASE TIME PP FLUID
	$(LINK) -shared -Wl,-soname,libgrains.so $(MPIFLAGS) $(ZFLAGS) -o $(LIBDIR)/libgrains.so $(OBJS_XML) $(OBJS_GEOM) $(OBJS_COMPO) $(OBJS_BASE) $(OBJS_TIME) $(OBJS_PP) $(OBJS_FLUID) 
	@echo 'Compilation is complete'
	@echo '--------------------------------'
	@echo 'Library libgrains.so generated !'

BASE   :  $(OBJS_BASE)
	@echo '*** Package Base compiled'
	@echo

COMPO  : $(OBJS_COMPO)
	@echo '*** Package Component compiled'
	@echo

FLUID  : $(OBJS_FLUID)
	@echo '*** Package Fluid compiled'
	@echo

GEOM   : $(OBJS_GEOM)
	@echo '*** Package Geometry compiled'
	@echo

PP     : $(OBJS_PP)
	@echo '*** Package PostProcessing compiled'
	@echo

TIME   : $(OBJS_TIME)
	@echo '*** Package TimeIntegration compiled'
	@echo
	
XML    : $(OBJS_XML)
	@echo '*** Package XML compiled'
	@echo	

MESS   :
	@echo '***********************************************************'
	@echo 'GRAINS compiled in' $(GRAINS_BITS_DEFAULT) 'bits with' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION) 'and' $(GRAINS_MPI_DISTRIB)-$(GRAINS_MPI_VERSION)
	@echo '***********************************************************'
	
START  :
	@echo 'MPI libraries = '$(LIBMPI_FOR_GRAINS___)
	@echo 'Compilation flags = $(GRAINS_GENCCFLAGS)$(GRAINS_MPICCFLAGS)$(COMPFLAGS)'
	@echo 'Macros = $(MACRO)'
	@echo
	@echo 'Start compiling'
	@echo	
	
CPINC  : 
	@echo 'Copy include files from subdirectories to' $(INCDIROUT)
	rm -f $(INCDIROUT)/*
	cp $(ALLINC) $(INCDIROUT)/	 	

$(OBJDIR)/%.o : %.cpp
	@echo 'Compiling '$<
	$(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@	

# DEPENDENCIES
UTLDIR  := $(GRAINS_HOME)/Tools/Perl
DEPEND  := $(INCDIR:%=-I%) $(GRAINS_SYSINC:%=-I%)

depend :
	@echo 'Generating dependencies'
	touch ./Maketmp
	makedepend -f ./Maketmp $(DEPEND) $(SRCFILES_WP_BASE) \
	$(SRCFILES_WP_COMPO) $(SRCFILES_WP_FLUID) $(SRCFILES_WP_GEOM) \
	$(SRCFILES_WP_PP) $(SRCFILES_WP_TIME) $(SRCFILES_WP_XML)
	$(UTLDIR)/makeadd.pl $(OBJDIR)/ ./Maketmp ./Makedepend
	rm -f ./Maketmp ./Maketmp.bak

include ./Makedepend      # DO NOT DELETE

# CLEAN
clean : 
	@echo 'Grains : Deleting object and library files'
	rm -f $(INCDIROUT)/*
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/*

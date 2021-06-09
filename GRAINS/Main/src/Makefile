# REPERTOIRES
EXEDIR = ../bin$(GRAINS_FULL_EXT)
OBJDIR = ../obj$(GRAINS_FULL_EXT)

# FICHIERS
MAINSRC = main.cpp
MAINOBJ = $(MAINSRC:%.cpp=$(OBJDIR)/%.o)



# DEPEND OF
GRAINS      = ../../Grains
GRAINSINC   = $(GRAINS)/inc -I$(OCTREE_HOME)/DLMFD
GRAINSLIB   = $(GRAINS)/lib$(GRAINS_FULL_EXT)
GRAINSOBJ   = $(GRAINS)/obj$(GRAINS_FULL_EXT)/*.o
GRAINSFLAGS = -L$(GRAINSLIB) -lgrains



# COMPILATEUR
# Note: '-D_GLIBCXX_USE_CXX11_ABI=0' added for compatibility with C++11 installed on ubuntu distribution.
ifeq ($(GRAINS_SERCOMPIL_ENV),GNU)
#  GCCFLAGS    = -g -W -Wno-long-long -Wno-ctor-dtor-privacy -Wno-unused-parameter
  GCCFLAGS    = -pedantic -W -Wno-long-long -Wno-ctor-dtor-privacy -Wno-unused-parameter -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11
else
  GCCFLAGS    =
endif

COMP  = $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) -m$(GRAINS_BITS_DEFAULT) $(GRAINS_COMPIL_OPT) -fPIC -Wall -Wwrite-strings -Wconversion -Wshadow -Wno-deprecated -Wno-comment $(GCCFLAGS)
LINK  = $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) -m$(GRAINS_BITS_DEFAULT) $(GRAINS_COMPIL_OPT) -fPIC -Wall -Wwrite-strings -Wconversion -Wshadow -Wno-deprecated -Wno-comment $(GCCFLAGS)

DEBUG  = -g

MACRO       = -D_XML -D_TIMER
XERCESINC   = $(XERCESC_ROOT)/include
XERCESLIB   = $(XERCESC_ROOT)/lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
XERCESFLAGS = -L$(XERCESLIB) -lxerces-c -lxerces-depdom

MPIINC	    = $(GRAINS_MPI_INCDIR)
MPILIB      = $(GRAINS_MPI_LIBDIR)
MPIFLAGS    = -L$(MPILIB) $(LIBMPI_FOR_GRAINS___)

 
INCFLAGS  = $(XERCESINC:%=-I%) $(GRAINSINC:%=-I%) $(MPIINC:%=-I%)
COMPFLAGS = $(DEBUG) $(INCFLAGS) -fPIC
LINKFLAGS = $(DEBUG)



# -----------------------------------------------------------------------------
# Ne pas modifier ci-dessous (G.FERRER)
.SILENT:

all : grains

# EXECUTABLE
grains : MESS $(EXEDIR)/grains
	@echo 'Logiciel grains construit !'

$(EXEDIR)/grains : app
	@echo 'Construction de grains'
	@echo 'Linkage command for generation of grains:' $(LINK) -o $@ $(MAINOBJ) $(LINKFLAGS) $(GRAINSFLAGS) $(MPIFLAGS) $(XERCESFLAGS)
	$(LINK) $(LINKFLAGS) -o $@ $(MAINOBJ) $(GRAINSFLAGS) $(MPIFLAGS) $(XERCESFLAGS)
	#original line executed: $(LINK)  $(LINKFLAGS) $(XERCESFLAGS) $(GRAINSFLAGS) $(MPIFLAGS)  -$(MAINOBJ) -o $@

MESS   :
	@echo '********************************************'
	@echo 'MAIN Compilation en' $(GRAINS_BITS_DEFAULT) 'bits avec' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
	@echo '********************************************'

# PACKAGE
app :  Grains $(MAINOBJ)
	@echo 'Generation de Grains'

Grains :
	echo $(COMP) $(COMPFLAGS) $(MACRO) $@
	( cd ../../Grains/src; make ; cd ../../Main/src )
	
$(OBJDIR)/%.o : %.cpp
	@echo 'Compilation de '$<
	$(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@
	echo $(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@
	
# DEPENDANCES
UTLDIR  = $(HOME)/Utils
INCLUDE = $(GRAINSINC)
DEPEND  = -I$(INCLUDE) -I$(OCTREE_HOME)/DLMFD $(SYSINC)

depend :
	@echo 'Generation des dependances'
	touch ./Maketmp
	makedepend -f ./Maketmp $(DEPEND) $(MAINSRC)
	$(UTLDIR)/makeadd.pl $(OBJDIR)/ ./Maketmp ./Makedepend
	rm -f ./Maketmp ./Maketmp.bak

include ./Makedepend      # DO NOT DELETE

# CLEAN
clean :
	( cd ../../Grains/src; make clean; cd ../../Main/src )
	@echo 'Main   : Suppression des fichiers .o'
	rm -f $(MAINOBJ)
	rm -f $(EXEDIR)/*

# HELP
help :
	@echo 'Usage  : gmake [target]'
	@echo 'Target : '
	@echo '  grains : Construction executable grains2D (defaut)'
	@echo '  depend : Generation des dependances '
	@echo '  clean  : Suppression des fichiers .o'

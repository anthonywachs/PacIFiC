# REPERTOIRES
EXEDIR = ../bin$(GRAINS_FULL_EXT)
OBJDIR = ../obj$(GRAINS_FULL_EXT)

# FICHIERS
MAINSRC = orientation.cpp
MAINOBJ = $(MAINSRC:%.cpp=$(OBJDIR)/%.o)

# DEPEND OF
GRAINS      = ../../../Grains
GRAINSINC   = $(GRAINS)/inc
GRAINSLIB   = $(GRAINS)/lib$(GRAINS_FULL_EXT)
GRAINSOBJ   = $(GRAINS)/obj$(GRAINS_FULL_EXT)/*.o
GRAINSFLAGS = -L$(GRAINSLIB) -lgrains

XERCESINC   = $(XERCESC_ROOT)/include
XERCESLIB   = $(XERCESC_ROOT)/lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
XERCESFLAGS = -L$(XERCESLIB) -lxerces-c -lxerces-depdom

# COMPILATEUR
COMP  = $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
LINK  = $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
DEBUG = -g

INCFLAGS  = $(XERCESINC:%=-I%) $(GRAINSINC:%=-I%)
COMPFLAGS = $(DEBUG) $(INCFLAGS) -fPIC
LINKFLAGS = $(DEBUG)


# -----------------------------------------------------------------------------
# Ne pas modifier ci-dessous (G.FERRER)
.SILENT:

all : pos

# EXECUTABLE
pos : MESS $(EXEDIR)/orientation
	@echo 'Logiciel orientation construit !'

$(EXEDIR)/orientation : $(MAINOBJ)
	@echo 'Construction de orientation'
	$(LINK) -o $@ $(MAINOBJ) $(LINKFLAGS) $(GRAINSFLAGS) $(XERCESFLAGS) 
	#original line executed was: $(LINK) $(LINKFLAGS) $(XERCESFLAGS) $(GRAINSFLAGS) $(MAINOBJ) -o $@

MESS   :
	@echo '*************************************'
	@echo 'Compilation en' $(GRAINS_BITS_DEFAULT) 'bits avec' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
	@echo '*************************************'

$(OBJDIR)/%.o : %.cpp
	@echo 'Compilation de '$<
	$(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@


# DEPENDANCES
UTLDIR  = $(HOME)/Utils
INCLUDE = 

DEPEND  = -I$(INCLUDE) $(SYSINC)

depend :
	@echo 'Generation des dependances'
	touch ./Maketmp
	makedepend -f ./Maketmp $(DEPEND) $(MAINSRC)
	$(UTLDIR)/makeadd.pl $(OBJDIR)/ ./Maketmp ./Makedepend
	rm -f ./Maketmp ./Maketmp.bak

include ./Makedepend      # DO NOT DELETE

# CLEAN
clean :
	rm -f $(MAINOBJ)
	rm -f $(EXEDIR)/*
  

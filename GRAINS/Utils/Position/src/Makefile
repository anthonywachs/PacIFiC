# REPERTOIRES
EXEDIR = ../bin$(GRAINS_FULL_EXT)
OBJDIR = ../obj$(GRAINS_FULL_EXT)

# FICHIERS
MAINSRC = position.cpp
MAINOBJ = $(MAINSRC:%.cpp=$(OBJDIR)/%.o)

# COMPILATEUR
COMP  = $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
LINK  = $(GRAINS_CPP) -m$(GRAINS_BITS_DEFAULT)
DEBUG = -g

INCFLAGS  = 
COMPFLAGS = $(DEBUG) $(INCFLAGS) -fPIC
LINKFLAGS = $(DEBUG)


# -----------------------------------------------------------------------------
# Ne pas modifier ci-dessous (G.FERRER)
.SILENT:

all : pos

# EXECUTABLE
pos : MESS $(EXEDIR)/position
	@echo 'Logiciel position construit !'

$(EXEDIR)/position : $(MAINOBJ)
	@echo 'Construction de position'
	$(LINK) $(LINKFLAGS) $(MAINOBJ) -o $@

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
  

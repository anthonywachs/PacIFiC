# REPERTOIRES
SRCDIR = .
INCDIR = ../inc
OBJDIR = ../obj$(GRAINS_FULL_EXT)
LIBDIR = ../lib$(GRAINS_FULL_EXT)

# FICHIERS
SRCS_APP = \
	Grains_BuilderFactory.cpp \
	Grains.cpp \
	GrainsCoupledWithFluid.cpp \
	GrainsCoupledWithFluidMPI.cpp \
	App.cpp \
	AppFluide_Puits.cpp \
	AppFluide_Drag.cpp \
	AppFluide_LubricationCorrection.cpp \
	AppFluide_Temperature.cpp \
	AppSec.cpp \
	InterfaceFluide_BuilderFactory.cpp \
	InterfaceFluide.cpp \
	InterfaceFluide2D.cpp \
	InterfaceFluide3D.cpp \
	EnsComposant.cpp \
	LinkedCell.cpp \
	AddedMass.cpp \
	Voisins.cpp \
	GrainsMPI.cpp \
	MPIWrapperGrains.cpp \
	MPIWrapperGrainsPer.cpp \
	computingtime.cpp \
	solvercomputingtime.cpp \
	PostProcessingWriter.cpp \
	Text_PostProcessingWriter.cpp \
	DLMFD_PostProcessingWriter.cpp \
	GMV_PostProcessingWriter.cpp \
	Paraview_PostProcessingWriter.cpp \
	CompFeatures_PostProcessingWriter.cpp \
	Matlab_PostProcessingWriter.cpp \
	SystemState.cpp \
	GrainsGeomTest.cpp \
	GrainsParameters.cpp \
	Grains_Exec.cpp \
	GrainsCompFeatures.cpp \
	GrainsPorosity.cpp \
	STL_PostProcessingWriter.cpp \
	PostProcessingWriter_BuilderFactory.cpp \
	GrainsMPITest.cpp
OBJS_APP = $(SRCS_APP:%.cpp=$(OBJDIR)/%.o)


SRCS_GRAINS = \
	Contact_BuilderFactory.cpp \
	Cellule.cpp \
	CineObstacle.cpp \
	CineObstacle_F.cpp \
	CineParticule.cpp \
	Cinematique_BuilderFactory.cpp \
	CohContact.cpp\
	CompObstacle.cpp \
	Composant.cpp \
	ContactLaw.cpp \
	CylindricalBoxObstacle.cpp \
	Forme.cpp \
	FormeVdW.cpp \
	LeapFrog_3D.cpp \
	LeapFrog_Sphere.cpp \
	MonObstacle.cpp \
	MyContact.cpp \
	ERContact.cpp \
	ERHContact.cpp \
	Obstacle_BuilderFactory.cpp \
	Obstacle.cpp \
	ObstacleAbsorbant.cpp \
	ObstaclePeriodique.cpp \
	ObstacleChargement.cpp \
	ObstacleChargement_F.cpp \
	Particule.cpp \
	ParticulePeriodique.cpp \
	CompParticule.cpp \
	CompParticulePeriodique.cpp \
	ElementParticule.cpp \
	PointContact.cpp \
	Torseur.cpp \
	Erreur.cpp \
	TimeIntegrator.cpp \
	SecondOrderLeapFrog.cpp \
	TimeIntegrator_BuilderFactory.cpp \
	FirstOrderExplicit.cpp \
	SecondOrderExplicit.cpp \
	SecondOrderAdamsBashforth.cpp \
	SecondOrderRungeKutta.cpp \
	ThirdOrderAdamsBashforth.cpp
OBJS_GRAINS = $(SRCS_GRAINS:%.cpp=$(OBJDIR)/%.o)


SRCS_2D = \
	Disque.cpp \
	LeapFrog_2D.cpp
OBJS_2D = $(SRCS_2D:%.cpp=$(OBJDIR)/%.o)


SRCS_SOLID = \
	Convex_BuilderFactory.cpp \
	BBox.cpp \
	Box.cpp \
	Cone.cpp \
	Convex.cpp \
	Cylinder.cpp \
	AntiCylinder.cpp \
	Plan.cpp \
	PointC.cpp \
	Polygon.cpp \
	Polyhedron.cpp \
	Polytope.cpp \
	Segment.cpp \
	Sphere.cpp \
	Transform.cpp \
	VertexBase.cpp \
	Plan_FiniC.cpp \
	Plan_InfiniC.cpp
OBJS_SOLID = $(SRCS_SOLID:%.cpp=$(OBJDIR)/%.o)


SRCS_3D = \
	Group3.cpp \
	IndexArray.cpp \
	Matrix.cpp \
	Point.cpp \
	Quaternion.cpp \
	Vecteur.cpp
OBJS_3D = $(SRCS_3D:%.cpp=$(OBJDIR)/%.o)

SRCS_UTL = ReaderXML.cpp WriterXML.cpp
OBJS_UTL = $(SRCS_UTL:%.cpp=$(OBJDIR)/%.o)

# COMPILATEUR
# Note: '-D_GLIBCXX_USE_CXX11_ABI=0' added (presumably) for compatibility ubuntu distribution
ifeq ($(GRAINS_SERCOMPIL_ENV),GNU)
  GCCFLAGS    = -pedantic -W -Wno-long-long -Wno-ctor-dtor-privacy -Wno-unused-parameter  -DMPICH_IGNORE_CXX_SEEK  -DMPICH_SKIP_MPICXX -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11
else
  GCCFLAGS    =
endif

COMP  = $(GRAINS_MPI_BINDIR)/$(GRAINS_MPICCC) -m$(GRAINS_BITS_DEFAULT) $(GRAINS_COMPIL_OPT) -fPIC -Wall -Wwrite-strings -Wconversion -Wshadow -Wno-deprecated -Wno-comment $(GCCFLAGS)
LINK  = $(GRAINS_MPI_BINDIR)/$(MPIBIN)/$(GRAINS_MPICCC) -m$(GRAINS_BITS_DEFAULT) $(GRAINS_COMPIL_OPT) -fPIC -Wall -Wwrite-strings -Wconversion -Wshadow -Wno-deprecated -Wno-comment $(GCCFLAGS)

DEBUG  = -g

MACRO       = -D_XML
XERCESINC   = $(XERCESC_ROOT)/include
XERCESLIB   = $(XERCESC_ROOT)/lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
XERCESFLAGS = -L$(XERCESLIB) -lxerces-c -lxerces-depdom

MPIINC	    = $(GRAINS_MPI_INCDIR)
MPILIB      = $(GRAINS_MPI_LIBDIR)
#MPIFLAGS    = -L$(MPILIB) $(LIBMPI_FOR_GRAINS___)

export LIBMPI_FOR_GRAINS___ = `echo ${GRAINS_MPI_LIBS} | sed 's%[^ ][^ ]*%-l&%g'`
MPIFLAGS    = -L$(MPILIB) $(LIBMPI_FOR_GRAINS___)

ZLIB = /usr/lib$(GRAINS_BITS_EXT)
ZFLAGS = -L$(ZLIB) -lz

INCFLAGS  = $(INCDIR:%=-I%) -I$(OCTREE_HOME)/DLMFD  $(XERCESINC:%=-I%) $(MPIINC:%=-I%)
COMPFLAGS = $(DEBUG) $(INCFLAGS) -fPIC

#------------------------------------------------------------------------------
# Ne pas modifier ci-dessous (G.FERRER)
.SILENT:

all : grains Basilisk
	@echo '-----------------------'

grains : MESS UTL 3D SOLID 2D GRAINS APP
	$(COMP) -o $(LIBDIR)/libgrains.so $(OBJS_APP) $(OBJS_GRAINS) $(OBJS_2D) $(OBJS_SOLID) $(OBJS_3D) $(OBJS_UTL) -shared -Wl,-soname,libgrains.so $(MPIFLAGS) $(ZFLAGS)
	#Original instruction executed was: $(COMP) -shared -Wl,-soname,libgrains.so $(MPIFLAGS) $(ZFLAGS) -o $(LIBDIR)/libgrains.so $(OBJS_APP) $(OBJS_GRAINS) $(OBJS_2D) $(OBJS_SOLID) $(OBJS_3D) $(OBJS_UTL)
	@echo 'Package GRAINS genere !'

Basilisk :
	$(COMP) $(COMPFLAGS) -shared $(SRCDIR)/InterfaceGrains.cpp -L$(LIBDIR) -lgrains $(XERCESFLAGS) -o $(LIBDIR)/libInterfaceGrains.so
	rm -f $(BASILISK_DLMFD)/InterfaceGrains.h
	ln -s $(GRAINS_ROOT)/inc/InterfaceGrains.h $(BASILISK_DLMFD)/InterfaceGrains.h

APP    :  $(OBJS_APP)
	@echo 'Package App genere'

GRAINS : $(OBJS_GRAINS)
	@echo 'Package Grains genere'

2D     : $(OBJS_2D)
	@echo 'Package specifique 2D genere'

SOLID  : $(OBJS_SOLID)
	@echo 'Package Solid genere'

3D     : $(OBJS_3D)
	@echo 'Package 3D genere'

UTL    : $(OBJS_UTL)
	@echo 'Package UTL genere'

MESS   :
	@echo '********************************************'
	@echo 'GRAINS Compilation en' $(GRAINS_BITS_DEFAULT) 'bits avec' $(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION)
	@echo '********************************************'
	@echo 'Librairies MPI = '$(LIBMPI_FOR_GRAINS___)

$(OBJDIR)/%.o : %.cpp
	@echo 'Compilation de '$<
	$(COMP) -o $@ -c $< $(COMPFLAGS) $(MACRO)
	#@echo $(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@
	#original instruction executed was: $(COMP) $(COMPFLAGS) $(MACRO) -c $< -o $@
	#Displays the paths where includes files are looked for:
	#echo | gcc -E $(COMPFLAGS) $(MACRO) -Wp,-v -

# DEPENDANCES
UTLDIR  = $(HOME)/Utils

DEPEND  = -I$(INCDIR) $(INCLUDE) $(SYSINC)

depend :
	@echo 'Generation des dependances'
	touch ./Maketmp
	makedepend -f ./Maketmp $(DEPEND) \
	$(SRCS_3D) $(SRCS_SOLID) $(SRCS_2D) $(SRCS_GRAINS) $(SRCS_APP)
	$(UTLDIR)/makeadd.pl $(OBJDIR)/ ./Maketmp ./Makedepend
	rm -f ./Maketmp ./Maketmp.bak

include ./Makedepend      # DO NOT DELETE

# CLEAN
clean :
	@echo 'Grains : Suppression des fichiers .o'
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/*
	rm -f -r $(OBJDIR)/ii_files

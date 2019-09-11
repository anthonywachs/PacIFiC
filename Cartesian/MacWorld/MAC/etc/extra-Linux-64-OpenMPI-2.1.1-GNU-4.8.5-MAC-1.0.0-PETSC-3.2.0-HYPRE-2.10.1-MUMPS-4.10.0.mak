# Common Makefile for Linux plateform runnning gcc
# Settings for Extra Components

ifeq ($(LINK_MAC),1)

WITH_MPI      = 1

##### PETSc external API
# if PETSc is enabled, uncomment the appropriate definition of PETSc_VERSION
WITH_PETSc    = 1
PETSc_VERSION = 3.2.0

WITH_MUMPS    = 0
WITH_INTEL    = 1


endif

WITH_ZLIB = 1
WITH_EXTENDED_MATH = 1


###############################################################################
ifeq ($(WITH_MUMPS),1)

ifeq ($(MAKE_MAC),1)
SRC += $(wildcard $(MAC_HOME)/ExternalAPI/MUMPS/src/*.cc)
CPPFLAGS += -I$(MAC_HOME)/ExternalAPI/MUMPS/include
endif

#PARMETIS = -L$(PETSC_DIR)/../MUMPS/parmetis-4.0.2/lib -lparmetis -L$(PETSC_DIR)/../MUMPS/parmetis-4.0.2/build/Linux-x86_64/libmetis -lmetis
#PARMETIS =

# The order of the linking with the static libraries for MUMPS matter
# Apparently it should be kept to : MUMPS - SCALAPACK - GFORTRAN

CPPFLAGS += -I$(MUMPS_DIR)/include 
LIBPATH  += $(MUMPS_DIR)/lib
LDLIBS   += -ldmumps-Linux-$(MACWORLD_FULL_EXT) -lmumps_common-Linux-$(MACWORLD_FULL_EXT) -lpord-Linux-$(MACWORLD_FULL_EXT)

LIBPATH  += $(MACWORLD_SCALAPACK_LIBDIR)
LDLIBS   += $(LIBSCALAPACK_FOR_MAC___)

LIBPATH  += $(MACWORLD_GFORTRAN_LIBDIR)
LDLIBS   += -l$(MACWORLD_GFORTRAN_LIBS)

WITH_MPI = 1
WITH_BLAS = 1
endif


###############################################################################
ifeq ($(WITH_PETSc),1)

ifeq ($(MAKE_MAC),1)
SRC += $(wildcard $(MAC_HOME)/ExternalAPI/PETSc_$(PETSc_VERSION)/src/*.cc)
CPPFLAGS += -I$(MAC_HOME)/ExternalAPI/PETSc_$(PETSc_VERSION)/include
endif

CPPFLAGS += -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIBPATH  += $(PETSC_DIR)/$(PETSC_ARCH)/lib
LDLIBS   += -lpetsc

# For Hypre
LIBPATH  += $(HYPRE_DIR)/$(HYPRE_ARCH)/lib
LDLIBS   += -lHYPRE

WITH_X11 = 1
WITH_BLAS = 1
endif


###############################################################################
ifeq ($(WITH_MPI),1)

# Mandatory definition MPIRUN
# this variable should contain the complete path to the mpirun command
# (without using any other variable)
MPIRUN = ${MACWORLD_MPI_BINDIR}/mpirun

ifeq ($(MAKE_MAC),1)
SRC += $(wildcard $(MAC_HOME)/ExternalAPI/MPI/src/*.cc)
CPPFLAGS += -I$(MAC_HOME)/ExternalAPI/MPI/include -DMPIRUN=\"$(MPIRUN)\"
endif

MPIPATH   = $(MACWORLD_MPI_ROOT)
CPPFLAGS += -I$(MACWORLD_MPI_INCDIR)
LIBPATH  += $(MACWORLD_MPI_LIBDIR)
LDLIBS   += $(LIBMPI_FOR_MAC___)

endif


###############################################################################
ifeq ($(WITH_X11),1)

X11PATH = /usr
CPPFLAGS += -I$(X11PATH)/include/X11
LIBPATH  += $(X11PATH)/lib$(MACWORLD_BITS_EXT)/
LDLIBS   += -lnsl -lXt -lX11 -lXmu
endif


###############################################################################
ifeq ($(WITH_SYSF77),1)
#LDLIBS   += -lg2c -lm
LDLIBS   += -lm
endif


###############################################################################
ifeq ($(WITH_EXTENDED_MATH),1)
CPPFLAGS += -DEXTENDED_MATH
endif


###############################################################################
ifeq ($(WITH_BLAS),1)
LIBPATH  += /usr/lib$(MACWORLD_BITS_EXT)
LDLIBS   += -lm
LIBPATH  += $(MACWORLD_BLAS_LIBDIR)
LDLIBS   += $(LIBLAS_FOR_MAC___)
LIBPATH  += $(MACWORLD_ATLAS_LIBDIR)
LDLIBS   += $(LIBATLAS_FOR_MAC___)
LIBPATH  += $(MACWORLD_LAPACK_LIBDIR)
LDLIBS   += $(LIBLAPACK_FOR_MAC___)
endif


###############################################################################
ifeq ($(WITH_ZLIB),1)
ZLIBPATH = /usr
LIBPATH  += $(ZLIBPATH)/lib$(MACWORLD_BITS_EXT)
CPPFLAGS += -I$(ZLIBPATH)/include -DZLIB
LDLIBS   += -lz
endif


###############################################################################
ifeq ($(WITH_INTEL),1)
LIBPATH  += $(MACWORLD_INTEL_LIBDIR)
LDLIBS   += $(LIBINTEL_FOR_MAC___)
endif

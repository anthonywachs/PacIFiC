# Common Makefile for Linux plateform runnning gcc
# Settings for Extra Components

ifeq ($(LINK_MAC),1)

WITH_MPI      = 1
WITH_PETSc    = 1
PETSc_VERSION = $(PETSC_VERSION)

endif

WITH_ZLIB = 1
WITH_EXTENDED_MATH = 1

ifeq ($(MACWORLD_SERCOMPIL_ENV),GNU)
  WITH_GNU = 1
else
  WITH_GNU = 0
endif

ifeq ($(MACWORLD_SERCOMPIL_ENV),Intel)
  WITH_INTEL = 1
else
  WITH_INTEL = 0
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
CPPFLAGS += -I$(MACWORLD_X11_INCDIR)
LIBPATH  += $(MACWORLD_X11_LIBDIR)
LDLIBS   += -lnsl -lXt -lX11 -lXmu

endif



###############################################################################
ifeq ($(WITH_EXTENDED_MATH),1)

CPPFLAGS += -DEXTENDED_MATH

endif



###############################################################################
ifeq ($(WITH_BLAS),1)

LIBPATH  += $(MACWORLD_M_LIBDIR)
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

LIBPATH  += $(MACWORLD_Z_LIBDIR)
CPPFLAGS += -I$(MACWORLD_Z_INCDIR) -DZLIB
LDLIBS   += -lz

endif



###############################################################################
ifeq ($(WITH_INTEL),1)

LIBPATH  += $(MACWORLD_INTEL_LIBDIR)
LDLIBS   += $(LIBINTEL_FOR_MAC___)

endif



###############################################################################
ifeq ($(WITH_GNU),1)

LIBPATH  += $(MACWORLD_GFORTRAN_LIBDIR)
LDLIBS   += $(LIBGNU_FOR_MAC___)

endif

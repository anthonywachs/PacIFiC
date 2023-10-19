# Common Makefile for Linux plateform runnning gcc
flagsdbg  = -g
flagsopt0 = ${MAC_OPT_FLAGS}
flagsopt1 = ${MAC_OPT_FLAGS}
flagsopt2 = ${MAC_OPT_FLAGS}

# setting the valid profile flags when needed
ifeq ($(WITH_PROFILE),1)
OPT += -pg
endif
# setting the valid test coverage flags when needed
ifeq ($(WITH_COVERAGE),1)
OPT += -fprofile-arcs -ftest-coverage
endif

CC       = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_C}
CXX      = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_CXX}
CPP      = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_C}
CXXFLAGS += $(OPT) -fPIC -Wall -Wno-long-long -Wno-ctor-dtor-privacy -pedantic -W -Wcast-qual -Wwrite-strings -Wconversion -Wshadow -Wno-unused-parameter
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -DOMPI_IGNORE_CXX_SEEK -DOMPI_SKIP_MPICXX
CXXFLAGS += -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11

FC       = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_F77}
FFLAGS  += $(OPT) -fPIC

LDFLAGS += $(OPT) -fPIC
LDFLAGS += -lm

LD.so       = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_CXX} -shared -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO = 

MKDEP.c  = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_C} -M $(CPPFLAGS)
MKDEP.cc = ${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_CXX} -M $(CPPFLAGS)

# Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -rpath -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) -Xlinker -rpath -Xlinker $(path))
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)

# Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))


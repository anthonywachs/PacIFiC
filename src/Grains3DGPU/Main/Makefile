#
# Environment
#

CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
NVCC = $(PACIFIC_NVCC)

MAIN_GPU_BINDIR = $(MAIN_GPU_BUILDDIR)/bin
MAIN_GPU_OBJDIR = $(MAIN_GPU_BUILDDIR)/obj

LIBGRAINS_GPU_INCDIR = $(LIBGRAINS_GPU_BUILDDIR)/include
LIBGRAINS_GPU_LIBDIR = $(LIBGRAINS_GPU_BUILDDIR)/lib

#
# NVCC + Compile + Link flags
#

include ../mk/nvcc_link_flags.mk

NVCC_ARCHFLAGS += -arch=sm_${PACIFIC_NVCC_ARCH}
NVCCFLAGS += -O3 -g -Xcompiler -fPIC -x cu -rdc=true --extended-lambda --expt-relaxed-constexpr

CPPFLAGS += -I$(LIBGRAINS_GPU_INCDIR)
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS)

LDFLAGS += -L$(LIBGRAINS_GPU_LIBDIR)

LDLIBS += -lGrains3DGPU
LDLIBS += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_CUDART_LFLAGS) $(PACIFIC_CURAND_LFLAGS)
NVCC_LDLIBS := $(call to_nvcc_ldflags,$(LDLIBS))

# 
# Target directories
#

LIBGRAINS_GPU_SO := $(LIBGRAINS_GPU_LIBDIR)/libGrains3DGPU.so
RPATH_FLAGS := -Xlinker -rpath -Xlinker $(LIBGRAINS_GPU_LIBDIR)

TARGET := $(MAIN_GPU_BINDIR)/grains3dgpu
SRC := src/main.cpp
OBJ := $(MAIN_GPU_OBJDIR)/main.o
DEP := $(OBJ:.o=.d)

#
# Production rules
#

.PHONY: all clean

all: $(TARGET)

$(MAIN_GPU_BINDIR):
	@mkdir -p $@

$(OBJ): $(SRC) | $(MAIN_GPU_BINDIR)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCC_ARCHFLAGS) $(CPPFLAGS) $(NVCCFLAGS) -MMD -MP -MF $(DEP) -MT $@ -c $< -o $@

$(TARGET): $(OBJ) $(LIBGRAINS_GPU_SO)
	$(NVCC) $(NVCC_ARCHFLAGS) $(RPATH_FLAGS) -o $@ $(OBJ) $(LDFLAGS) $(NVCC_LDLIBS)
	@echo -e '-----------------------\n''Software grains3dgpu built !\n''-----------------------'

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)

-include $(DEP)

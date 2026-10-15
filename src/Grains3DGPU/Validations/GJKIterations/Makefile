#
# Environment (typically provided by parent Makefiles)
#

VALIDATIONS_GPU_BUILDDIR ?= $(PACIFIC_BUILDDIR_ABS)/Grains3DGPU/Validations
LIBGRAINS_GPU_BUILDDIR ?= $(PACIFIC_BUILDDIR_ABS)/Grains3DGPU
LIBGRAINS_GPU_SRCDIR ?= $(PACIFIC_SRCDIR_ABS)/Grains3DGPU/Grains

#
# Derived environment
#

VALIDATIONS_GPU_BINDIR := $(VALIDATIONS_GPU_BUILDDIR)/bin

LIBGRAINS_GPU_OBJDIR := $(LIBGRAINS_GPU_BUILDDIR)/obj
LIBGRAINS_GPU_INCDIR := $(LIBGRAINS_GPU_BUILDDIR)/include
LIBGRAINS_GPU_LIBDIR := $(LIBGRAINS_GPU_BUILDDIR)/lib
OBJLIB_MANIFEST := $(LIBGRAINS_GPU_LIBDIR)/libGrains3DGPU_objects.mk

GJK_ITERATIONS_BINDIR := $(VALIDATIONS_GPU_BINDIR)
GJK_ITERATIONS_OBJDIR := $(VALIDATIONS_GPU_BUILDDIR)/GJKIterations/obj

GJK_ITERATIONS_DATADIR := $(VALIDATIONS_GPU_BUILDDIR)/GJKIterations/data
 

#
# Build flags
#

NVCC = $(PACIFIC_NVCC)
include ../../mk/nvcc_link_flags.mk

NVCC_ARCHFLAGS += -arch=sm_${PACIFIC_NVCC_ARCH}
NVCCFLAGS += -O3 -g -Xcompiler -fPIC -x cu -rdc=true --extended-lambda --expt-relaxed-constexpr

CPPFLAGS += -I$(LIBGRAINS_GPU_INCDIR)
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS) $(PACIFIC_CURAND_CFLAGS)
CPPFLAGS += -D_XML

LDLIBS += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_CUDART_LFLAGS) $(PACIFIC_CURAND_LFLAGS)
NVCC_LDLIBS := $(call to_nvcc_ldflags,$(LDLIBS))

SRC := main.cpp
OBJ := $(GJK_ITERATIONS_OBJDIR)/main.o
DEP := $(OBJ:.o=.d)
TARGET := $(GJK_ITERATIONS_BINDIR)/GJKIterationsTest

#
# Production rules
#

.PHONY: all test clean clean-data clean-all objlib-manifest

all: objlib-manifest $(TARGET)

objlib-manifest:
	@$(MAKE) -C $(LIBGRAINS_GPU_SRCDIR) objlib \
		LIBGRAINS_GPU_BUILDDIR=$(LIBGRAINS_GPU_BUILDDIR)

$(GJK_ITERATIONS_BINDIR):
	@mkdir -p $@

$(GJK_ITERATIONS_OBJDIR):
	@mkdir -p $@

$(GJK_ITERATIONS_DATADIR):
	@mkdir -p $@

$(OBJ): $(SRC) | $(GJK_ITERATIONS_OBJDIR)
	$(NVCC) $(NVCC_ARCHFLAGS) $(CPPFLAGS) $(NVCCFLAGS) -MMD -MP -MF $(DEP) -MT $@ -c $< -o $@

$(TARGET): $(OBJ) objlib-manifest | $(GJK_ITERATIONS_BINDIR) $(GJK_ITERATIONS_DATADIR)
	@objs="$$(sed -n 's/^LIBGRAINS_GPU_OBJECTS :=//p' $(OBJLIB_MANIFEST))"; \
	if [ -z "$$objs" ]; then \
		echo "Error: could not load object list from $(OBJLIB_MANIFEST)"; \
		exit 1; \
	fi; \
	$(NVCC) $(NVCC_ARCHFLAGS) -o $@ $(OBJ) $$objs $(NVCC_LDLIBS)
	@echo -e '--------------------------------\n''Validation GJKIterationsTest built !\n''--------------------------------'

test: $(TARGET)
	$(TARGET) $(TEST_ARGS)

clean:
	rm -rf $(GJK_ITERATIONS_OBJDIR) $(TARGET)

clean-data:
	rm -rf $(GJK_ITERATIONS_DATADIR)

clean-all: clean clean-data

-include $(DEP)

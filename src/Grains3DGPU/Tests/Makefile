# Compilation & linking flags

THIS_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
LIBGRAINS_GPU_SRCDIR ?= $(abspath $(THIS_DIR)/../Grains)

LIBGRAINS_GPU_INCDIR = $(LIBGRAINS_GPU_BUILDDIR)/include
LIBGRAINS_GPU_LIBDIR = $(LIBGRAINS_GPU_BUILDDIR)/lib
LIBGRAINS_GPU_OBJDIR = $(LIBGRAINS_GPU_BUILDDIR)/obj
OBJLIB_MANIFEST = $(LIBGRAINS_GPU_LIBDIR)/libGrains3DGPU_objects.mk

TEST_GPU_OBJDIR = $(TESTS_GPU_BUILDDIR)/obj
TEST_GPU_BINDIR = $(TESTS_GPU_BUILDDIR)/bin

NVCC = $(PACIFIC_NVCC)
include ../mk/nvcc_link_flags.mk

NVCC_ARCHFLAGS += -arch=sm_${PACIFIC_NVCC_ARCH}
NVCCFLAGS += -O3 -g -Xcompiler -fPIC -x cu -rdc=true --extended-lambda --expt-relaxed-constexpr

CPPFLAGS += -I$(LIBGRAINS_GPU_INCDIR)
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS) $(PACIFIC_CURAND_CFLAGS)
CPPFLAGS += $(PACIFIC_GTEST_CFLAGS) $(PACIFIC_GTEST_MAIN_CFLAGS)
CPPFLAGS += -DTEST_BUILD=1

LDLIBS += $(PACIFIC_GTEST_MAIN_LFLAGS) $(PACIFIC_GTEST_LFLAGS)
LDLIBS += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_CUDART_LFLAGS) $(PACIFIC_CURAND_LFLAGS)
NVCC_LDLIBS := $(call to_nvcc_ldflags,$(LDLIBS))

# Target directories
TEST_GPU_BINDIR ?= $(GRAINS_GPU_BINDIR)
TEST_GPU_OBJDIR ?= $(GRAINS_GPU_OBJDIR)/tests

TARGET := $(TEST_GPU_BINDIR)/grains3dgpu_tests

# Sources
TEST_SRCS := \
	main.cpp \
	geometry/test_vector3.cpp \
	geometry/test_matrix3.cpp \
	geometry/test_quaternion.cpp \
	geometry/RotationMathTest.cpp \
	collision/test_gjk.cpp \
	collision/test_obb.cpp \
	collision/test_obc.cpp \
	collision/test_collision_detection.cpp \
	collision/test_neighborlist.cpp \
	contactForce/test_contact_table_host.cpp \
	contactForce/test_contact_table_device.cpp

TEST_SRCS := $(foreach f,$(TEST_SRCS),$(if $(wildcard $(f)),$(f),))
TEST_OBJS := $(patsubst %.cpp,$(TEST_GPU_OBJDIR)/%.o,$(TEST_SRCS))
TEST_DEPS := $(TEST_OBJS:.o=.d)

# -----------------------------------------------------------------------------
# Production rules
.PHONY: all run clean objlib-manifest

all: $(TARGET)

objlib-manifest:
	@$(MAKE) -C $(LIBGRAINS_GPU_SRCDIR) objlib \
		LIBGRAINS_GPU_BUILDDIR=$(LIBGRAINS_GPU_BUILDDIR)

$(TEST_GPU_BINDIR):
	@mkdir -p $@

$(TEST_GPU_OBJDIR):
	@mkdir -p $@

$(TEST_GPU_OBJDIR)/%.o: %.cpp | $(TEST_GPU_OBJDIR)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCC_ARCHFLAGS) $(CPPFLAGS) $(NVCCFLAGS) -MMD -MP -MF $(@:.o=.d) -MT $@ -c $< -o $@

$(TARGET): $(TEST_OBJS) objlib-manifest | $(TEST_GPU_BINDIR)
	@objs="$$(sed -n 's/^LIBGRAINS_GPU_OBJECTS :=//p' $(OBJLIB_MANIFEST))"; \
	if [ -z "$$objs" ]; then \
		echo "Error: could not load object list from $(OBJLIB_MANIFEST)"; \
		exit 1; \
	fi; \
	$(NVCC) $(NVCC_ARCHFLAGS) -o $@ $(TEST_OBJS) $$objs -L$(LIBGRAINS_GPU_LIBDIR) $(NVCC_LDLIBS)
	@echo -e '-----------------------\n''Tests grains3dgpu_tests built !\n''-----------------------'

run: $(TARGET)
	$(TARGET)

clean:
	rm -rf $(TEST_GPU_OBJDIR) $(TARGET)

-include $(TEST_DEPS)

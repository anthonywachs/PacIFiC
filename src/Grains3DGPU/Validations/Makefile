#
# Environment (overwritten by command line args)
#

VALIDATIONS_GPU_BUILDDIR := $(PACIFIC_BUILDDIR_ABS)/Grains3DGPU/Validations
LIBGRAINS_GPU_BUILDDIR := $(PACIFIC_BUILDDIR_ABS)/Grains3DGPU

#
# Derived Environment
#

VALIDATIONS_GPU_BINDIR = $(VALIDATIONS_GPU_BUILDDIR)/bin

LIBGRAINS_GPU_OBJDIR = $(LIBGRAINS_GPU_BUILDDIR)/obj
LIBGRAINS_GPU_INCDIR = $(LIBGRAINS_GPU_BUILDDIR)/include
LIBGRAINS_GPU_LIBDIR = $(LIBGRAINS_GPU_BUILDDIR)/lib

#
# Production rules
#

all: collision-detection gjk-iterations
	@echo "All validation tools built successfully!"

#
# Collision Detection Module
#

COLLISION_DETECTION_MODULE_MAKE = @$(MAKE) -C CollisionDetectionModule \
	VALIDATIONS_GPU_BUILDDIR=$(VALIDATIONS_GPU_BUILDDIR) \
	LIBGRAINS_GPU_BUILDDIR=$(LIBGRAINS_GPU_BUILDDIR)

collision-detection-module:
	@echo "Building collision detection validation tools..."
	@$(COLLISION_DETECTION_MODULE_MAKE) all

#
# GJK Iteration 
#

GJK_ITERATIONS_MAKE = @$(MAKE) -C GJKIterations \
	VALIDATIONS_GPU_BUILDDIR=$(VALIDATIONS_GPU_BUILDDIR) \
	LIBGRAINS_GPU_BUILDDIR=$(LIBGRAINS_GPU_BUILDDIR)

gjk-iterations:
	@echo "Building GJK iterations validation tools..."
	@$(GJK_ITERATIONS_MAKE) all

# test-gjk-iterations: gjk-iterations
# 	@echo "Running GJK iterations validation test..."
# 	@$(GJK_ITERATIONS_MAKE) test

#
# Other
#

clean:
	@echo "Cleaning collision detection validation tools..."
	@$(COLLISION_DETECTION_MODULE_MAKE) clean
	@echo "Cleaning GJK iterations validation tools..."
	@$(GJK_ITERATIONS_MAKE) clean
	@echo "All validation tools cleaned!"

help:
	@echo "GrainsGPU Validations Makefile"
	@echo "=============================="
	@echo ""
	@echo "Available targets:"
	@echo "  all                      - Build all validation tools"
	@echo "  collision-detection      - Build collision detection validation tools"
	@echo "  gjk-iterations           - Build gjk iteration validation tools"
	@echo "  clean                    - Clean all validation tools"
	@echo "  help                     - Show this help message"
	@echo ""
	@echo "Individual validation directories:"
	@echo "  CollisionDetection/ - GJK performance comparison tools"
	@echo "  GJKIterations/      - GJK iteration count vs distance study"

# SUBMODULES := $(shell git submodule status --recursive | awk '{print $$2}')

# Programs
GIT := git
MKDIR := mkdir
CMAKE := cmake 
CP := cp
LN := ln

MAKEFLAGS += --no-print-directory

PACIFIC_MPI_SKIP_CPPFLAGS ?= -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1

PACIFIC_THIRDPARTY_MAKE = @$(MAKE) -C third_party \
  PACIFIC_ROOT_ABS="$(PACIFIC_ROOT_ABS)" \
  PACIFIC_THIRDPARTY_BUILDDIR_ABS=${PACIFIC_THIRDPARTY_BUILDDIR_ABS} \
  PACIFIC_THIRDPARTY_INSTALLDIR_ABS=${PACIFIC_THIRDPARTY_INSTALLDIR_ABS} \
  MKDIR="$(MKDIR)" CP="$(CP)" LN="$(LN)" GIT="$(GIT)" CMAKE="$(CMAKE)"  

PACIFIC_SRC_MAKE = @$(MAKE) -C $(PACIFIC_SRCDIR_ABS) \
  PACIFIC_ROOT_ABS="$(PACIFIC_ROOT_ABS)" \
  PACIFIC_SRCDIR_ABS="$(PACIFIC_SRCDIR_ABS)" \
  PACIFIC_BUILDDIR_ABS="$(PACIFIC_BUILDDIR_ABS)" \
  PACIFIC_MPI_SKIP_CPPFLAGS="$(PACIFIC_MPI_SKIP_CPPFLAGS)" \
  MKDIR="$(MKDIR)" CP="$(CP)" LN="$(LN)" GIT="$(GIT)" CMAKE="$(CMAKE)"

#
# Check environment for thirdparty deps requested
#

SRC_THIRDPARTY_DEPS := 
GRAINS_THIRDPARTY_DEPS :=
GRAINS_GPU_THIRDPARTY_DEPS :=
MAC_THIRDPARTY_DEPS :=
FLUID_THIRDPARTY_DEPS :=
OCTREE_THIRDPARTY_DEPS :=

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_BASILISK_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-basilisk
OCTREE_THIRDPARTY_DEPS += third_party-basilisk
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_HDF5_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-hdf5
OCTREE_THIRDPARTY_DEPS += third_party-hdf5
MAC_THIRDPARTY_DEPS += third_party-hdf5
MAC_THIRDPARTY_DEPS += third_party-hdf5
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_XERCESC_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-xercesc
GRAINS_THIRDPARTY_DEPS += third_party-xercesc
GRAINS_GPU_THIRDPARTY_DEPS += third_party-xercesc
FLUID_THIRDPARTY_DEPS += third_party-xercesc
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_ZLIB_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-zlib
GRAINS_THIRDPARTY_DEPS += third_party-zlib
GRAINS_GPU_THIRDPARTY_DEPS += third_party-zlib
MAC_THIRDPARTY_DEPS += third_party-zlib
FLUID_THIRDPARTY_DEPS += third_party-zlib
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_PETSC_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-petsc
MAC_THIRDPARTY_DEPS += third_party-petsc
FLUID_THIRDPARTY_DEPS += third_party-petsc
endif

#
# Main targets
#

.PHONY: all docs docs-develop

all: grains3d mac fluid octree
all: grains3d mac fluid octree

# src: builddir $(SRC_THIRDPARTY_DEPS)
# 	$(PACIFIC_SRC_MAKE) all

grains3d: builddir $(GRAINS_THIRDPARTY_DEPS)
grains3d: builddir $(GRAINS_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) grains

grains3dgpu: builddir venv $(GRAINS_GPU_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) grainsgpu

testsgpu: builddir $(GRAINS_GPU_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) testsgpu

mac: builddir $(MAC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) mac

fluid: builddir $(FLUID_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) fluid

octree: $(OCTREE_THIRDPARTY_DEPS)

docs:
	@$(MAKE) -C "$(PACIFIC_ROOT_ABS)/docs" build

docs-develop:
	@$(MAKE) -C "$(PACIFIC_ROOT_ABS)/docs" develop

#
# Git Submodules
#
 
submodule: submodule-basilisk submodule-hdf5 submodule-zlib submodule-xercesc submodule-petsc submodule-gtest
submodule-clean: submodule-basilisk-clean submodule-hdf5-clean submodule-zlib-clean submodule-xercesc-clean submodule-petsc-clean submodule-gtest-clean

submodule-basilisk: 
	@$(GIT) submodule update --init --checkout third_party/basilisk/basilisk_submodule
submodule-basilisk-clean:
	@$(GIT) submodule deinit -f -- third_party/basilisk/basilisk_submodule 

submodule-hdf5: 
	@$(GIT) submodule update --init third_party/hdf5/hdf5_submodule
submodule-hdf5-clean:
	@$(GIT) submodule deinit -f -- third_party/hdf5/hdf5_submodule 

submodule-zlib: 
	@$(GIT) submodule update --init third_party/zlib/zlib_submodule
submodule-zlib-clean:
	@$(GIT) submodule deinit -f -- third_party/zlib/zlib_submodule 

submodule-xercesc: 
	@$(GIT) submodule update --init third_party/xercesc/xercesc_submodule
submodule-xercesc-clean:
	@$(GIT) submodule deinit -f -- third_party/xercesc/xercesc_submodule 

submodule-petsc: 
	@$(GIT) submodule update --init third_party/petsc/petsc_submodule
submodule-petsc-clean:
	@$(GIT) submodule deinit -f -- third_party/petsc/petsc_submodule 

submodule-gtest: 
	@$(GIT) submodule update --init third_party/gtest/gtest_submodule
submodule-gtest-clean:
	@$(GIT) submodule deinit -f -- third_party/gtest/gtest_submodule 

#
# thirdparty
#

builddir:
	$(MKDIR) -p $(PACIFIC_BUILDDIR_ABS)

third_party: builddir submodule
	$(PACIFIC_THIRDPARTY_MAKE) all

third_party-basilisk: builddir submodule-basilisk
	$(PACIFIC_THIRDPARTY_MAKE) third_party-basilisk

third_party-hdf5: builddir submodule-hdf5
	$(PACIFIC_THIRDPARTY_MAKE) third_party-hdf5

third_party-zlib: builddir submodule-zlib
	$(PACIFIC_THIRDPARTY_MAKE) third_party-zlib

third_party-xercesc: builddir submodule-xercesc
	$(PACIFIC_THIRDPARTY_MAKE) third_party-xercesc

third_party-petsc: builddir submodule-petsc
	$(PACIFIC_THIRDPARTY_MAKE) third_party-petsc

third_party-gtest: builddir submodule-gtest
	$(PACIFIC_THIRDPARTY_MAKE) third_party-gtest

clean: submodule-clean
	$(PACIFIC_SRC_MAKE) clean
	$(PACIFIC_THIRDPARTY_MAKE) clean

#
# Python
#

PACIFIC_VENV_ABS ?= $(PACIFIC_BUILDDIR_ABS)/.venv
PACIFIC_VENV_PIP ?= $(PACIFIC_VENV_ABS)/bin/pip

venv: builddir
	python3 -m venv ${PACIFIC_VENV_ABS}
	$(PACIFIC_VENV_PIP) install --upgrade pip
	$(PACIFIC_VENV_PIP) install -r ${PACIFIC_ROOT_ABS}/requirements.txt

#
# Help
#

RESET  := \033[0m
BOLD   := \033[1m
BLUE   := \033[34m
GREEN  := \033[32m
CYAN   := \033[36m
YELLOW := \033[33m

help:
	@printf "Targets:\n"
	@printf "\t$(GREEN)make all   $(RESET)\n"
	@printf "\t$(GREEN)make clean $(RESET)\n"
	@printf "\t$(GREEN)make grains3d$(RESET)\n"
	@printf "\t$(GREEN)make grains3dgpu$(RESET)\n" 
	@printf "\t$(GREEN)make mac   $(RESET)\n"
	@printf "\t$(GREEN)make fluid $(RESET)\n"
	@printf "\t$(CYAN)make docs        $(RESET)\n"
	@printf "\t$(CYAN)make docs-develop$(RESET)\n"
	@printf "\t$(RED)make third_party         $(RESET)\n"
	@printf "\t$(RED)make third_party-basilisk$(RESET)\n"
	@printf "\t$(RED)make third_party-hdf5    $(RESET)\n"
	@printf "\t$(RED)make third_party-zlib    $(RESET)\n"
	@printf "\t$(RED)make third_party-xercesc $(RESET)\n"
	@printf "\t$(RED)make third_party-gtest   $(RESET)\n"
	@printf "\t$(YELLOW)make submodule$(RESET)\n"
	@printf "\t$(YELLOW)make submodule-clean$(RESET)\n"

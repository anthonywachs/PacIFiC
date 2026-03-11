# SUBMODULES := $(shell git submodule status --recursive | awk '{print $$2}')

# Programs
GIT := git
MKDIR := mkdir
CMAKE := cmake 
CP := cp
LN := ln
CC := gcc 

PACIFIC_THIRDPARTY_MAKE = $(MAKE) -C third_party \
  PACIFIC_ROOT_ABS="$(PACIFIC_ROOT_ABS)" \
  PACIFIC_THIRDPARTY_BUILDDIR_ABS=${PACIFIC_THIRDPARTY_BUILDDIR_ABS} \
  PACIFIC_THIRDPARTY_INSTALLDIR_ABS=${PACIFIC_THIRDPARTY_INSTALLDIR_ABS} \
  MKDIR="$(MKDIR)" CP="$(CP)" LN="$(LN)" GIT="$(GIT)" CMAKE="$(CMAKE)"  

PACIFIC_SRC_MAKE = $(MAKE) -C $(PACIFIC_SRCDIR_ABS) \
  PACIFIC_ROOT_ABS="$(PACIFIC_ROOT_ABS)" \
  PACIFIC_SRCDIR_ABS="$(PACIFIC_SRCDIR_ABS)" \
  PACIFIC_BUILDDIR_ABS="$(PACIFIC_BUILDDIR_ABS)" \
  MKDIR="$(MKDIR)" CP="$(CP)" LN="$(LN)" GIT="$(GIT)" CMAKE="$(CMAKE)"

#
# Check environment for thirdparty deps requested
#

SRC_THIRDPARTY_DEPS :=

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_BASILISK_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-basilisk
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_HDF5_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-hdf5
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_XERCESC_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-xercesc
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_ZLIB_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-zlib
endif

ifneq ($(filter 1 true TRUE yes YES on ON,$(PACIFIC_PETSC_USE_THIRDPARTY)),)
SRC_THIRDPARTY_DEPS += third_party-petsc
endif

.PHONY: all docs docs-develop

all: src

src: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) all

grains: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) grains

mac: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) mac

fluid: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) fluid

docs:
	$(MAKE) -C "$(PACIFIC_ROOT_ABS)/docs" build

docs-develop:
	$(MAKE) -C "$(PACIFIC_ROOT_ABS)/docs" develop

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
	@printf "\t$(GREEN)make grains$(RESET)\n"
	@printf "\t$(GREEN)make mac   $(RESET)\n"
	@printf "\t$(GREEN)make fluid $(RESET)\n"
	@printf "\t$(CYAN)make docs        $(RESET)\n"
	@printf "\t$(CYAN)make docs-develop$(RESET)\n"
	@printf "\t$(RED)make third_party         $(RESET)\n"
	@printf "\t$(RED)make third_party-basilisk$(RESET)\n"
	@printf "\t$(RED)make third_party-hdf5    $(RESET)\n"
	@printf "\t$(RED)make third_party-zlib    $(RESET)\n"
	@printf "\t$(RED)make third_party-xercesc $(RESET)\n"
	@printf "\t$(YELLOW)make submodules$(RESET)\n"

submodules:
	@$(GIT) submodule update --init --recursive
	@$(GIT) submodule update --init --recursive --checkout third_party/basilisk/basilisk_submodule


# submodules:

builddir:
	$(MKDIR) -p $(PACIFIC_BUILDDIR_ABS)

third_party: builddir submodules
	$(PACIFIC_THIRDPARTY_MAKE) all

third_party-basilisk: builddir submodules
	$(PACIFIC_THIRDPARTY_MAKE) third_party-basilisk

third_party-hdf5: builddir submodules
	$(PACIFIC_THIRDPARTY_MAKE) third_party-hdf5

third_party-zlib: builddir submodules
	$(PACIFIC_THIRDPARTY_MAKE) third_party-zlib

third_party-xercesc: builddir submodules
	$(PACIFIC_THIRDPARTY_MAKE) third_party-xercesc

third_party-petsc: builddir submodules
	$(PACIFIC_THIRDPARTY_MAKE) third_party-petsc

clean:
	$(PACIFIC_SRC_MAKE) clean
	$(PACIFIC_THIRDPARTY_MAKE) clean

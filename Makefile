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

.PHONY: all docs docs-develop

all: src

src: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) all

grains: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) grains

mac: builddir $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) mac

fluid: q $(SRC_THIRDPARTY_DEPS)
	$(PACIFIC_SRC_MAKE) fluid

docs:
	$(MAKE) -C "$(PACIFIC_ROOT_ABS)/docs" build

docs-develop:
	$(MAKE) -C "$(PACIFIC_ROOT_ABS)/docs" develop

# help:
# 	@echo "Targets:"
# 	@echo "\tmake submodules"

# submodules:
# 	@$(GIT) submodule sync --recursive 
# 	@$(GIT) submodule update --init --recursive

submodules:

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

clean:
	$(PACIFIC_SRC_MAKE) clean
	$(PACIFIC_THIRDPARTY_MAKE) clean

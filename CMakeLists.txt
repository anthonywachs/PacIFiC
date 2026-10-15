# ================================================================================
#
# Top level project settings
#
# ================================================================================

cmake_minimum_required(VERSION 3.12.0...4.0.0)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(defaults)

# ================================================================================
#
# Set the project version
#
# ================================================================================

find_package(Git QUIET)

set(PACIFIC_FALLBACK_VERSION "0.0.1")

# User-facing overrides (set via -D...=... or cmake-gui)
set(PACIFIC_VERSION_FULL_OVERRIDE   "" CACHE STRING "Override full version string (highest priority).")
set(PACIFIC_VERSION_SEMVER_OVERRIDE "" CACHE STRING "Override SemVer (e.g. 1.2.3). Used if FULL override is empty.")

# Defaults
set(PACIFIC_VERSION_SEMVER "${PACIFIC_FALLBACK_VERSION}")
set(PACIFIC_VERSION_FULL   "${PACIFIC_FALLBACK_VERSION}")

# 1) User-provided FULL override (highest priority)
if(PACIFIC_VERSION_FULL_OVERRIDE)
  set(PACIFIC_VERSION_FULL "${PACIFIC_VERSION_FULL_OVERRIDE}")

  # SemVer: prefer explicit semver override, else try to parse from FULL, else fallback
  if(PACIFIC_VERSION_SEMVER_OVERRIDE)
    set(PACIFIC_VERSION_SEMVER "${PACIFIC_VERSION_SEMVER_OVERRIDE}")
  else()
    string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" _m "${PACIFIC_VERSION_FULL}")
    if(CMAKE_MATCH_1)
      set(PACIFIC_VERSION_SEMVER "${CMAKE_MATCH_1}")
    endif()
  endif()

# 2) User-provided SemVer override (next priority)
elseif(PACIFIC_VERSION_SEMVER_OVERRIDE)
  set(PACIFIC_VERSION_SEMVER "${PACIFIC_VERSION_SEMVER_OVERRIDE}")
  set(PACIFIC_VERSION_FULL   "${PACIFIC_VERSION_SEMVER_OVERRIDE}")

# 3) Git-provided (if enabled and available)
elseif(PACIFIC_VERSION_USE_GIT AND GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --tags --match "v[0-9]*" --abbrev=12 --dirty --always
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  if(GIT_DESCRIBE)
    string(REGEX REPLACE "^v" "" PACIFIC_VERSION_FULL "${GIT_DESCRIBE}")
    string(REGEX MATCH "^([0-9]+\\.[0-9]+\\.[0-9]+)" _m "${PACIFIC_VERSION_FULL}")
    if(CMAKE_MATCH_1)
      set(PACIFIC_VERSION_SEMVER "${CMAKE_MATCH_1}")
    else()
      # No reachable SemVer tag; keep fallback semver but still set a meaningful FULL
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      if(GIT_SHA)
        set(PACIFIC_VERSION_FULL "${PACIFIC_VERSION_SEMVER}-g${GIT_SHA}")
      endif()
    endif()
  endif()
# Fallback already set above
endif()

# ================================================================================
#
# Project definition
#
# ================================================================================

message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "PacIFiC version ${PACIFIC_VERSION_FULL}")

project(PacIFiC 
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES C CXX
)

# Avoid deprecated MPI C++ bindings in mpi.h across all targets.
add_compile_definitions(
  MPICH_SKIP_MPICXX=1
  OMPI_SKIP_MPICXX=1
)

# ==========================================================
#
# Custom Build Configurations
#
# ==========================================================

# ==========================================================
# Custom build type: RelWithWarnings
#
# This sets some flags similar to those used in packaging 
# for various Linux distributions.
# ==========================================================

set(CMAKE_CXX_FLAGS_RELWITHWARNINGS "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELWITHWARNINGS  "${CMAKE_C_FLAGS_RELEASE}")

set(_WARN_FLAGS
  "-Werror=return-type \
   -Werror=format-security \
   -Werror=format=2 \
   -Werror=nonnull \
   -Werror=array-bounds \
   -Werror=stringop-overflow \
   -Werror=stringop-truncation \
   -Werror=alloc-size-larger-than=9223372036854775807 \
   -Werror=address \
   -Werror=maybe-uninitialized \
   -Werror=mismatched-new-delete \
   -Werror=overloaded-virtual \
   -Wunused-but-set-variable \
   -Wunused-variable \
   -Wparentheses \
   -Wcomment \
   -Wsign-compare \
   -Wunused-function \
   -Wdangling-else \
   -Wmisleading-indentation \
   -Wreorder"
)

set(_SANITIZER_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer")

string(APPEND CMAKE_CXX_FLAGS_RELWITHWARNINGS
       " ${_WARN_FLAGS} ${_SANITIZER_FLAGS}")
string(APPEND CMAKE_C_FLAGS_RELWITHWARNINGS
       " ${_WARN_FLAGS} ${_SANITIZER_FLAGS}")

mark_as_advanced(CMAKE_CXX_FLAGS_RELWITHWARNINGS CMAKE_C_FLAGS_RELWITHWARNINGS)


# ------------------------------------------------------------------------------
# Handle optional MPI install prefix: This is entirely used for Fedora/RHEL 
# packaging, where MPI libraries are installed in a non-standard location
# ------------------------------------------------------------------------------

include(GNUInstallDirs)

if(DEFINED CMAKE_MPI_INSTALL_LIBDIR AND CMAKE_MPI_INSTALL_LIBDIR)
  set(PACIFIC_MPI_INSTALL_LIBDIR "${CMAKE_MPI_INSTALL_LIBDIR}")
else()
  set(PACIFIC_MPI_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
endif()


if(DEFINED CMAKE_MPI_INSTALL_BINDIR AND CMAKE_MPI_INSTALL_BINDIR)
  set(PACIFIC_MPI_INSTALL_BINDIR "${CMAKE_MPI_INSTALL_BINDIR}")
else()
  set(PACIFIC_MPI_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
endif()

if(DEFINED CMAKE_MPI_INSTALL_INCLUDEDIR AND CMAKE_MPI_INSTALL_INCLUDEDIR)
  set(PACIFIC_MPI_INSTALL_INCLUDEDIR "${CMAKE_MPI_INSTALL_INCLUDEDIR}")
else()
  set(PACIFIC_MPI_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
endif()


# ================================================================================
#
# Build and install third-party dependencies as part of CMake build
#
# By default, these are OFF. They can be enabled by setting the corresponding
# PACIFIC_BUILD_THIRD_PARTY_<LIBRARY> option to ON in ccmake or using 
# -DPACIFIC_BUILD_THIRD_PARTY_<LIBRARY>=ON
# 
# ================================================================================

add_subdirectory(third_party)  

if(PACIFIC_BUILD_TESTS)
  enable_testing()
endif()

# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

add_subdirectory(src)

# ================================================================================
#
# Docs
#
# ================================================================================

if(PACIFIC_BUILD_DOCS)
  add_subdirectory(docs)
endif()

# ================================================================================
#
# Environment File
#
# ================================================================================

include(PacificEnv)
 

if (PACIFIC_DEV_ENV) 
  if (OCTREE_BASILISK_PROVIDER STREQUAL "SYSTEM")
    pacific_generate_dev_env_script(
      OUTPUT  "${CMAKE_BINARY_DIR}/pacific-env.sh"
      TARGETS Grains3D_main Grains3D_position Grains3D_rotationmatrix Grains3D_shapefile fluid
      # EXTRA_LIB_DIRS "${CMAKE_BINARY_DIR}/src/Grains3D/Grains"
    )
  else()
    pacific_generate_dev_env_script(
      OUTPUT  "${CMAKE_BINARY_DIR}/pacific-env.sh"
      TARGETS Grains3D_main Grains3D_position Grains3D_rotationmatrix Grains3D_shapefile fluid qcc
      # EXTRA_LIB_DIRS "${CMAKE_BINARY_DIR}/src/Grains3D/Grains"
    )
  endif()
  message(STATUS "Environment Script: ${CMAKE_BINARY_DIR}/pacific-env.sh")
endif()


 

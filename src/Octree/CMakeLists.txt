cmake_minimum_required(VERSION 3.20)
project(Octree LANGUAGES C)

# ------------------------------------------------------------------------------
# Core deps
# ------------------------------------------------------------------------------

if (PACIFIC_USE_MPI)
  find_package(MPI REQUIRED)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_HDF5)
else()
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 REQUIRED COMPONENTS C HL)
endif()

include(GNUInstallDirs)
include(FetchContent)

# ------------------------------------------------------------------------------
# Basilisk / qcc provider selection
#
#   SYSTEM   : use qcc found on PATH
#   VENDORED : use user-provided Basilisk source dir via 
#              -DBASILISK_DIR=/path/to/basilisk, or git submodule if 
#              -DUSE_SUBMODULES=ON
#   GIT      : fetch Basilisk from git using CMake's FetchContent
# ------------------------------------------------------------------------------
 
set(OCTREE_BASILISK_PROVIDER "GIT" CACHE STRING "How to provide Basilisk/qcc: SYSTEM, VENDORED, GIT")
set_property(CACHE OCTREE_BASILISK_PROVIDER PROPERTY STRINGS SYSTEM VENDORED GIT)

# What the wrapper should execute. Default is plain 'qcc' (works if we install qcc to bindir).
set(QCC_PROGRAM "qcc")

if(OCTREE_BASILISK_PROVIDER STREQUAL "SYSTEM")

  find_program(QCC_PROGRAM qcc REQUIRED)
  message(STATUS "Using system qcc: ${QCC_PROGRAM}")

elseif(OCTREE_BASILISK_PROVIDER STREQUAL "VENDORED")

  # Define the cache variable ONCE
  set(BASILISK_DIR "" CACHE PATH
      "Path to a Basilisk source tree (must contain CMakeLists.txt)")

  # Compute effective directory
  set(_basilisk_dir "${BASILISK_DIR}")

  if(PACIFIC_USE_SUBMODULES AND _basilisk_dir STREQUAL "")
    set(_basilisk_dir "${CMAKE_SOURCE_DIR}/third_party/basilisk/basilisk_submodule")
    message(STATUS "PACIFIC_USE_SUBMODULES=ON: using Basilisk submodule at ${_basilisk_dir}")
  endif()

  if(_basilisk_dir STREQUAL "")
    message(FATAL_ERROR
      "OCTREE_BASILISK_PROVIDER=VENDORED requires:\n"
      "  -DBASILISK_DIR=/path/to/basilisk\n"
      "or\n"
      "  -DUSE_SUBMODULES=ON (expects submodule checkout)")
  endif()

  if(NOT EXISTS "${_basilisk_dir}/CMakeLists.txt")
    message(FATAL_ERROR
      "Basilisk dir '${_basilisk_dir}' is missing CMakeLists.txt.\n"
      "If using submodules, run:\n"
      "  git submodule update --init --recursive")
  endif()

  # Correct FetchContent override
  set(FETCHCONTENT_SOURCE_DIR_BASILISK "${_basilisk_dir}" CACHE PATH "" FORCE)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK ON CACHE BOOL "" FORCE)
  mark_as_advanced(
    FETCHCONTENT_SOURCE_DIR_BASILISK
    FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK
  )
  FetchContent_Declare(Basilisk)
  FetchContent_MakeAvailable(Basilisk)

elseif(OCTREE_BASILISK_PROVIDER STREQUAL "GIT")
  FetchContent_Declare(
    Basilisk
    GIT_REPOSITORY https://github.com/tortotubus/pacific-basilisk.git
  )
  FetchContent_MakeAvailable(Basilisk)

  mark_as_advanced(
    FETCHCONTENT_SOURCE_DIR_BASILISK
    FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK
  )
else()
  message(FATAL_ERROR
    "Unknown OCTREE_BASILISK_PROVIDER='${OCTREE_BASILISK_PROVIDER}'. Use SYSTEM, VENDORED, or GIT.")
endif()



# ------------------------------------------------------------------------------
# Octree headers
# ------------------------------------------------------------------------------
set(OCTREE_DLMFD_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/DLMFD"        CACHE PATH "Octree DLMFD directory")
set(OCTREE_VTKHYPERTREE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/VTKHyperTree" CACHE PATH "Octree VTKHyperTree directory")
set(OCTREE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/Octree")
mark_as_advanced(OCTREE_DLMFD_DIR OCTREE_VTKHYPERTREE_DIR OCTREE_INSTALL_CMAKEDIR)

# Install

include(CMakePackageConfigHelpers)

add_library(Octree_headers INTERFACE)
add_library(PacIFiC::Octree ALIAS Octree_headers)

target_include_directories(Octree_headers INTERFACE
  $<BUILD_INTERFACE:${OCTREE_DLMFD_DIR}>
  $<BUILD_INTERFACE:${OCTREE_VTKHYPERTREE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Octree>
)

install(
  DIRECTORY "${OCTREE_DLMFD_DIR}/"
  COMPONENT Octree_Devel
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Octree/DLMFD"
  FILES_MATCHING PATTERN "*.h"
)

install(
  DIRECTORY "${OCTREE_VTKHYPERTREE_DIR}/"
  COMPONENT Octree_Devel
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Octree/VTKHyperTree"
  FILES_MATCHING PATTERN "*.h"
)

install(TARGETS Octree_headers
  EXPORT OctreeTargets
  COMPONENT Octree_Devel
)

install(EXPORT OctreeTargets
  FILE OctreeTargets.cmake
  NAMESPACE PacIFiC::
  DESTINATION "${OCTREE_INSTALL_CMAKEDIR}"
  COMPONENT Octree_Devel
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/OctreeConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/OctreeConfig.cmake"
  INSTALL_DESTINATION "${OCTREE_INSTALL_CMAKEDIR}"
)

# write_basic_package_version_file(
#   "${CMAKE_CURRENT_BINARY_DIR}/OctreeConfigVersion.cmake"
#   VERSION "${PROJECT_VERSION}"
#   COMPATIBILITY SameMajorVersion
# )

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/OctreeConfig.cmake"
  # "${CMAKE_CURRENT_BINARY_DIR}/OctreeConfigVersion.cmake"
  DESTINATION "${OCTREE_INSTALL_CMAKEDIR}"
  COMPONENT Octree_Devel
)

# --- Install helper CMake modules -------------------------------------------------
install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/octree_add_executable.cmake"
  DESTINATION "${OCTREE_INSTALL_CMAKEDIR}"
  COMPONENT Octree_Devel
)



# ------------------------------------------------------------------------------
# qcc wrapper script install
# ------------------------------------------------------------------------------
# For SYSTEM, QCC_PROGRAM is an absolute path (from find_program).
# For VENDORED/GIT, we install 'qcc' into our bindir, so calling "qcc" is fine.
# set(QCC_INSTALL_LOCATION "${QCC_PROGRAM}")

# set(QCC_EXTRA_OPTS
#   "-lm "
#   "-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/Octree/DLMFD "
#   "-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/Octree/VTKHyperTree "
#   "-L${CMAKE_INSTALL_FULL_LIBDIR} "
#   "-Wl,-rpath,${CMAKE_INSTALL_FULL_LIBDIR}"
#   "-lGrains ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES}"
# )

# configure_file(
#   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qcc-pacific.in"
#   "${CMAKE_CURRENT_BINARY_DIR}/qcc-pacific"
#   @ONLY
# )

# install(
#   PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/qcc-pacific"
#   DESTINATION "${CMAKE_INSTALL_BINDIR}"
# )

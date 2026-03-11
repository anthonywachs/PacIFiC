
cmake_minimum_required(VERSION 3.12.0...4.0.0)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(FLUID
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES CXX
)

#
# Collect FLUID sources
#

set(FLUID_SOURCES "")
set(FLUID_INCLUDE_DIRS "")

#
# FluidSolidInterface
#

list(APPEND FLUID_INCLUDE_DIRS $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FluidSolidInterface/include>)
file(GLOB FLUID_FSI_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/FluidSolidInterface/src/*.cc"
)
list(SORT FLUID_FSI_SOURCES)
list(APPEND FLUID_SOURCES ${FLUID_FSI_SOURCES})

#
# Misc
#

list(APPEND FLUID_INCLUDE_DIRS $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Misc/include>)
file(GLOB FLUID_MISC_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Misc/src/*.cc"
)
list(SORT FLUID_MISC_SOURCES)
list(APPEND FLUID_SOURCES ${FLUID_MISC_SOURCES})

#
# PostProcessing
#

list(APPEND FLUID_INCLUDE_DIRS $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include>)
file(GLOB FLUID_POSTPROCESSING_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/src/*.cc"
)
list(SORT FLUID_POSTPROCESSING_SOURCES)
list(APPEND FLUID_SOURCES ${FLUID_POSTPROCESSING_SOURCES})

#
# Solvers/DirectionSplitting
#

list(APPEND FLUID_INCLUDE_DIRS $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DirectionSplitting/include>)
file(GLOB FLUID_DS_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DirectionSplitting/src/*.cc"
)
list(SORT FLUID_DS_SOURCES)
list(APPEND FLUID_SOURCES ${FLUID_DS_SOURCES})

#
# Solvers/DLMFD
#

list(APPEND FLUID_INCLUDE_DIRS $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DLMFD/include>)
file(GLOB FLUID_DLMFD_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DLMFD/src/*.cc"
)
list(SORT FLUID_DLMFD_SOURCES)
list(APPEND FLUID_SOURCES ${FLUID_DLMFD_SOURCES})

#
# Solvers/Viscoplastic_HeatTransfer
#

list(APPEND FLUID_INCLUDE_DIRS $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Viscoplastic_HeatTransfer/include>)
file(GLOB FLUID_VPH_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Viscoplastic_HeatTransfer/src/*.cc"
)
list(SORT FLUID_VPH_SOURCES)
list(APPEND FLUID_SOURCES ${FLUID_VPH_SOURCES})

#
# Build FLUID binary
#
 
add_executable(fluid ${FLUID_SOURCES})
 

target_include_directories(fluid
  PUBLIC
    ${FLUID_INCLUDE_DIRS}
)

target_link_libraries(fluid 
  PUBLIC
    PacIFiC::MAC
    PacIFiC::Grains
)


# ================================================================================
#
# Install 
#
# ================================================================================

include(GNUInstallDirs)
# include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS fluid
  # EXPORT FLUIDTargets
  RUNTIME DESTINATION ${PACIFIC_MPI_INSTALL_BINDIR}
  COMPONENT FLUID_Runtime
)

# ---- 4. Generate CMake package configuration --------------------------------------
# set(FLUID_INSTALL_CMAKEDIR ${PACIFIC_MPI_INSTALL_LIBDIR}/cmake/FLUID)

# configure_package_config_file(
#   ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FLUIDConfig.cmake.in
#   ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfig.cmake
#   INSTALL_DESTINATION ${FLUID_INSTALL_CMAKEDIR}
# )

# write_basic_package_version_file(
#   ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfigVersion.cmake
#   VERSION ${PROJECT_VERSION}
#   COMPATIBILITY SameMajorVersion
# )

# ---- 5. Install the package + targets ---------------------------------------------
# install(EXPORT FLUIDTargets
#   FILE FLUIDTargets.cmake
#   NAMESPACE PacIFiC::
#   DESTINATION ${FLUID_INSTALL_CMAKEDIR}
#   COMPONENT FLUID_Devel
# )

# install(FILES
#   ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfig.cmake
#   ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfigVersion.cmake
#   DESTINATION ${FLUID_INSTALL_CMAKEDIR}
#   COMPONENT FLUID_Devel
# )

# ---- 6. Export from the build tree for development -------------------------------
# export(EXPORT FLUIDTargets
#   FILE ${CMAKE_CURRENT_BINARY_DIR}/FLUIDTargets.cmake
#   NAMESPACE PacIFiC::
# )

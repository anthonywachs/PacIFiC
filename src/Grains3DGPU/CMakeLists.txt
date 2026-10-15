cmake_minimum_required(VERSION 3.12.0...4.0.0)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(Grains3DGPU
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES CUDA CXX 
)


# ================================================================================
#
# USER PROVIDED DEPENDENCIES
#
# ================================================================================

if (PACIFIC_USE_MPI)
  find_package(MPI REQUIRED COMPONENTS C)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_ZLIB)
else()
  find_package(ZLIB REQUIRED)
endif()

if (PACIFIC_BUILD_THIRD_PARTY_XERCESC)
else()
  find_package(XercesC REQUIRED)
endif()

if (PACIFIC_USE_CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

if (PACIFIC_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if (PACIFIC_USE_GTEST)
  find_package(GTest REQUIRED)
endif()


# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

add_subdirectory(Grains)
add_subdirectory(Main) 
add_subdirectory(Validations) 


if(PACIFIC_BUILD_TESTS)
  set(GRAINS3D_GPU_SOURCE_ROOT_DIR "${PROJECT_SOURCE_DIR}")
  set(GRAINS3D_GPU_BINARY_ROOT_DIR "${PROJECT_BINARY_DIR}")
  enable_testing()
  add_subdirectory(Tests)
endif()


# ================================================================================
#
# Install 
#
# ================================================================================

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS 
  Grains3DGPU 
  Grains3DGPU_main
  EXPORT Grains3DGPUTargets
  LIBRARY DESTINATION ${PACIFIC_MPI_INSTALL_LIBDIR}      # libGrains3DGPU.so
  COMPONENT Grains3DGPU_Runtime
  RUNTIME DESTINATION ${PACIFIC_MPI_INSTALL_BINDIR}      # grains3dgpu
  COMPONENT Grains3DGPU_Runtime
  INCLUDES DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR} # for consumers
)

# ---- 2. Install all public headers -------------------------------------------------
install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Base/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/CollisionDetection/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Component/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/ContactForce/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Geometry/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/PostProcessing/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/TimeIntegration/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/XML3/include/
  DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}/Grains3DGPU
  COMPONENT Grains3DGPU_Devel
  FILES_MATCHING PATTERN "*.h*" # match .h, .hh, .hpp
)

# ---- 3. Install data/config files --------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Grains/Dtd/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/Grains3DGPU/Dtd
  COMPONENT Grains3D_Data
  FILES_MATCHING PATTERN "*.dtd"
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Grains/Tools/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/Grains3DGPU/Tools
  COMPONENT Grains3D_Data
  FILES_MATCHING PATTERN "*.exec"
)

# ---- 4. Generate CMake package configuration --------------------------------------
# This makes the project discoverable via find_package(Grains3D CONFIG)
set(GRAINS3D_GPU_INSTALL_CMAKEDIR ${PACIFIC_MPI_INSTALL_LIBDIR}/cmake/Grains3DGPU)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Grains3DGPUConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/Grains3DGPUConfig.cmake
  INSTALL_DESTINATION "${GRAINS3D_GPU_INSTALL_CMAKEDIR}"
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/Grains3DGPUConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# ---- 5. Install the package + targets ---------------------------------------------
install(EXPORT Grains3DGPUTargets
  FILE Grains3DGPUTargets.cmake
  NAMESPACE PacIFiC::
  DESTINATION ${GRAINS3D_GPU_INSTALL_CMAKEDIR}
  COMPONENT Grains3DGPU_Devel
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/Grains3DGPUConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/Grains3DGPUConfigVersion.cmake
  DESTINATION ${GRAINS3D_GPU_INSTALL_CMAKEDIR}
  COMPONENT Grains3DGPU_Devel
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT Grains3DGPUTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/Grains3DGPUTargets.cmake
  NAMESPACE PacIFiC::
)

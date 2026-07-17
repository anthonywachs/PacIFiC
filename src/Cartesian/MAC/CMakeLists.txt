cmake_minimum_required(VERSION 3.12.0...4.0.0) 

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(MAC
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES CXX
)

# ================================================================================
#
# USER PROVIDED DEPENDENCIES
#
# ================================================================================

find_package(MPI REQUIRED COMPONENTS C)

if (BUILD_THIRD_PARTY_ZLIB)
  # find_package(ZLIB REQUIRED)
else()
  find_package(ZLIB REQUIRED)
endif()

if (BUILD_THIRD_PARTY_HDF5)
else()
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 COMPONENTS C HL REQUIRED)
endif()

# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================
 
set(MAC_SOURCES "")
set(MAC_INCLUDE_DIRS "")

# ArrayVector
list(APPEND MAC_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ArrayVector/include>")
file(GLOB MAC_ARRAYVECTOR_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/ArrayVector/src/*.cc"
)
list(SORT MAC_ARRAYVECTOR_SOURCES)
list(APPEND MAC_SOURCES
  ${MAC_ARRAYVECTOR_SOURCES}
)

# Base
list(APPEND MAC_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Base/include>")
file(GLOB MAC_BASE_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Base/src/*.cc"
)
list(SORT MAC_BASE_SOURCES)
list(APPEND MAC_SOURCES
  ${MAC_BASE_SOURCES}
)

# ExternalAPI/MPI
list(APPEND MAC_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ExternalAPI/MPI/include>")
file(GLOB MAC_MPI_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/ExternalAPI/MPI/src/*.cc"
)
list(SORT MAC_MPI_SOURCES)
list(APPEND MAC_SOURCES
  ${MAC_MPI_SOURCES}
)

# ExternalAPI/PETSc_3.24.3
list(APPEND MAC_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ExternalAPI/PETSc_3.24.3/include>")
file(GLOB MAC_PETSC_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/ExternalAPI/PETSc_3.24.3/src/*.cc"
)
list(SORT MAC_PETSC_SOURCES)
list(APPEND MAC_SOURCES
  ${MAC_PETSC_SOURCES}
)

# FiniteVolume
list(APPEND MAC_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolume/include>")
file(GLOB MAC_FINITEVOLUME_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolume/src/*.cc"
)
list(SORT MAC_FINITEVOLUME_SOURCES)
list(APPEND MAC_SOURCES
  ${MAC_FINITEVOLUME_SOURCES}
)


# LinearAlgebra
list(APPEND MAC_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/include>")
file(GLOB MAC_LA_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/src/*.cc"
)
list(SORT MAC_LA_SOURCES)
list(APPEND MAC_SOURCES
  ${MAC_LA_SOURCES}
)

add_library(MAC SHARED ${MAC_SOURCES})
add_library(PacIFiC::MAC ALIAS MAC)

if(BUILD_THIRD_PARTY_PETSC)
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(PETSC REQUIRED PETSc)

  target_include_directories(MAC PRIVATE ${PETSC_INCLUDE_DIRS})
  target_link_directories(MAC PRIVATE ${PETSC_LIBRARY_DIRS})
  target_link_libraries(MAC PRIVATE ${PETSC_LIBRARIES})
  target_compile_options(MAC PRIVATE ${PETSC_CFLAGS_OTHER})
endif()


set_target_properties(MAC
  PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

set(MAC_LEVEL 0 CACHE STRING "Compile definition value for LEVEL in MAC")
set_property(CACHE MAC_LEVEL PROPERTY STRINGS "0" "1" "2")

target_compile_definitions(MAC
  PRIVATE
    LEVEL=${MAC_LEVEL}
    MPIRUN="${MPIEXEC_EXECUTABLE}"
)

message(STATUS "mpirun/mpiexec : ${MPIEXEC_EXECUTABLE}")

target_link_libraries(MAC 
  PUBLIC
    MPI::MPI_C
    ZLIB::ZLIB
    hdf5::hdf5
    hdf5::hdf5_hl
)

include(GNUInstallDirs)

target_include_directories(MAC
  PUBLIC
    ${MAC_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/MAC>
)


# ================================================================================
#
# Install 
#
# ================================================================================

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS MAC
  EXPORT MACTargets
  LIBRARY DESTINATION ${PACIFIC_MPI_INSTALL_LIBDIR}     
  COMPONENT MAC_Runtime
  RUNTIME DESTINATION ${PACIFIC_MPI_INSTALL_BINDIR}      
  COMPONENT MAC_Runtime
  INCLUDES DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}  
)

# ---- 2. Install all public headers -------------------------------------------------
# If your project has multiple include roots, add them all.
install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/ArrayVector/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Base/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/ExternalAPI/MPI/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/FiniteVolume/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/include/
  DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}/MAC
  COMPONENT MAC_Devel
  FILES_MATCHING PATTERN "*.h*"      # match .h, .hh, .hpp
)

# ---- 4. Generate CMake package configuration --------------------------------------
set(MAC_INSTALL_CMAKEDIR ${PACIFIC_MPI_INSTALL_LIBDIR}/cmake/MAC)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MACConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfig.cmake
  INSTALL_DESTINATION ${MAC_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# ---- 5. Install the package + targets ---------------------------------------------
install(EXPORT MACTargets
  FILE MACTargets.cmake
  NAMESPACE PacIFiC::
  DESTINATION ${MAC_INSTALL_CMAKEDIR}
  COMPONENT MAC_Devel
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MACConfigVersion.cmake
  DESTINATION ${MAC_INSTALL_CMAKEDIR}
  COMPONENT MAC_Devel
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT MACTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/MACTargets.cmake
  NAMESPACE PacIFiC::
)

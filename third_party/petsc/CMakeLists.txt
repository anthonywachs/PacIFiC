include(ExternalProject)
include(GNUInstallDirs)

set(PETSC_VERSION "3.25.3" CACHE STRING "PETSc version to fetch when USE_SUBMODULES=OFF")
set(PETSC_LIBDIR "lib" CACHE STRING "PETSc install library directory relative to PETSC_PREFIX")
set(PETSC_CONFIGURE_OPTIONS
  "--prefix=<INSTALL_DIR>"
  "--with-shared-libraries=1"
  "--with-debugging=0"
  "--with-fortran-bindings=0"
  CACHE STRING "Arguments passed to PETSc configure"
)

if(USE_SUBMODULES)
  set(_petsc_source_dir "${CMAKE_SOURCE_DIR}/third_party/petsc/petsc_submodule")
  set(_petsc_download_args "")
else()
  set(_petsc_source_dir "${CMAKE_BINARY_DIR}/_deps/petsc-src")
  set(_petsc_download_args
    GIT_REPOSITORY https://gitlab.com/petsc/petsc.git
    GIT_TAG        "v${PETSC_VERSION}"
  )
endif()

set(_petsc_lib "${PETSC_PREFIX}/${PETSC_LIBDIR}/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}")

ExternalProject_Add(petsc_ep
  ${_petsc_download_args}

  SOURCE_DIR  "${_petsc_source_dir}"
  BINARY_DIR  "${_petsc_source_dir}"
  INSTALL_DIR "${PETSC_PREFIX}"

  CONFIGURE_COMMAND
    ${CMAKE_COMMAND} -E env
      PETSC_DIR=<SOURCE_DIR>
      <SOURCE_DIR>/configure
      ${PETSC_CONFIGURE_OPTIONS}

  BUILD_COMMAND
    ${CMAKE_COMMAND} -E env
      PETSC_DIR=<SOURCE_DIR>
      ${CMAKE_MAKE_PROGRAM} -C <SOURCE_DIR> all

  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E env
      PETSC_DIR=<SOURCE_DIR>
      ${CMAKE_MAKE_PROGRAM} -C <SOURCE_DIR> install

  UPDATE_COMMAND ""
  UPDATE_DISCONNECTED 1
  BUILD_ALWAYS OFF
  BUILD_BYPRODUCTS
    "${_petsc_lib}"
)

file(MAKE_DIRECTORY "${PETSC_PREFIX}/include")
file(MAKE_DIRECTORY "${PETSC_PREFIX}/${PETSC_LIBDIR}")

add_library(PETSc::PETSc UNKNOWN IMPORTED GLOBAL)
set_target_properties(PETSc::PETSc PROPERTIES
  IMPORTED_LOCATION "${_petsc_lib}"
  INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${PETSC_PREFIX}/include>"
)

add_dependencies(PETSc::PETSc petsc_ep)

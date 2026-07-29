include(ExternalProject)
include(GNUInstallDirs)

set(PETSC_VERSION "3.25.3" CACHE STRING "PETSc version to fetch when USE_SUBMODULES=OFF")
set(PACIFIC_PETSC_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "PETSc install prefix")
set(PACIFIC_PETSC_LIBDIR "lib" CACHE STRING "PETSc install library directory relative to PACIFIC_PETSC_PREFIX")
set(PETSC_CONFIGURE_OPTIONS
  "--prefix=<INSTALL_DIR>"
  "--with-shared-libraries=1"
  "--with-debugging=0"
  "--with-fortran-bindings=0"
  CACHE STRING "Arguments passed to PETSc configure"
)

function(pacific_petsc_detect_source_version source_dir output_variable)
  set(_petsc_version_header "${source_dir}/include/petscversion.h")

  if(NOT EXISTS "${_petsc_version_header}")
    message(FATAL_ERROR
      "Could not detect PETSc version: '${_petsc_version_header}' does not exist.")
  endif()

  foreach(_component MAJOR MINOR SUBMINOR)
    set(_line_var "_petsc_${_component}_line")
    file(STRINGS "${_petsc_version_header}" ${_line_var}
      REGEX "^#define[ \t]+PETSC_VERSION_${_component}[ \t]+[0-9]+")

    if(NOT ${_line_var})
      message(FATAL_ERROR
        "Could not detect PETSc ${_component} version from '${_petsc_version_header}'.")
    endif()

    string(REGEX REPLACE
      "^#define[ \t]+PETSC_VERSION_${_component}[ \t]+([0-9]+).*$"
      "\\1"
      _value
      "${${_line_var}}")
    set(_petsc_${_component} "${_value}")
  endforeach()

  set(${output_variable}
    "${_petsc_MAJOR}.${_petsc_MINOR}.${_petsc_SUBMINOR}"
    PARENT_SCOPE)
endfunction()

if(USE_SUBMODULES)
  set(_petsc_source_dir "${CMAKE_SOURCE_DIR}/third_party/petsc/petsc_submodule")
  set(_petsc_download_args "")

  pacific_petsc_detect_source_version("${_petsc_source_dir}" _petsc_source_version)
  set(PETSC_VERSION "${_petsc_source_version}" CACHE STRING
    "PETSc version to fetch when USE_SUBMODULES=OFF" FORCE)
  message(STATUS
    "Using PETSc submodule at ${_petsc_source_dir} (version ${PETSC_VERSION})")
else()
  set(_petsc_source_dir "${CMAKE_BINARY_DIR}/_deps/petsc-src")
  set(_petsc_download_args
    GIT_REPOSITORY https://gitlab.com/petsc/petsc.git
    GIT_TAG        "v${PETSC_VERSION}"
  )
endif()

set(_petsc_lib "${PACIFIC_PETSC_PREFIX}/${PACIFIC_PETSC_LIBDIR}/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}")

ExternalProject_Add(petsc_ep
  ${_petsc_download_args}

  SOURCE_DIR  "${_petsc_source_dir}"
  BINARY_DIR  "${_petsc_source_dir}"
  INSTALL_DIR "${PACIFIC_PETSC_PREFIX}"
  STAMP_DIR   "${CMAKE_CURRENT_BINARY_DIR}/petsc_ep-stamp-${PETSC_VERSION}"

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

file(MAKE_DIRECTORY "${PACIFIC_PETSC_PREFIX}/include")
file(MAKE_DIRECTORY "${PACIFIC_PETSC_PREFIX}/${PACIFIC_PETSC_LIBDIR}")

add_library(PETSc::PETSc UNKNOWN IMPORTED GLOBAL)
set_target_properties(PETSc::PETSc PROPERTIES
  IMPORTED_LOCATION "${_petsc_lib}"
  INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${PACIFIC_PETSC_PREFIX}/include>"
)

add_dependencies(PETSc::PETSc petsc_ep)

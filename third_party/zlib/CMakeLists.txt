include(ExternalProject)
include(GNUInstallDirs)

if (PACIFIC_USE_SUBMODULES)
  ExternalProject_Add(zlibng_ep
    UPDATE_COMMAND ""
    SOURCE_DIR     "${CMAKE_SOURCE_DIR}/third_party/zlib/zlib_submodule"
    BINARY_DIR     "${CMAKE_BINARY_DIR}/_deps/zlib-ng-build"
    INSTALL_DIR    "${ZLIBNG_PREFIX}"
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DZLIB_COMPAT=ON
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
    BUILD_BYPRODUCTS
      "${ZLIBNG_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libz.so"
      "${ZLIBNG_PREFIX}/include/zlib.h"
  )
else()
  ExternalProject_Add(zlibng_ep
    GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
    GIT_TAG        2.2.4
    SOURCE_DIR     "${CMAKE_BINARY_DIR}/_deps/zlib-ng-src"
    BINARY_DIR     "${CMAKE_BINARY_DIR}/_deps/zlib-ng-build"
    INSTALL_DIR    "${ZLIBNG_PREFIX}"
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DZLIB_COMPAT=ON
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
    BUILD_BYPRODUCTS
      "${ZLIBNG_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libz.so"
      "${ZLIBNG_PREFIX}/include/zlib.h"
  )
endif()

file(MAKE_DIRECTORY "${ZLIBNG_PREFIX}/include")
file(MAKE_DIRECTORY "${ZLIBNG_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

add_library(ZLIB::ZLIB UNKNOWN IMPORTED GLOBAL)
set(_zlib_lib "${ZLIBNG_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libz.so")
set_target_properties(ZLIB::ZLIB PROPERTIES
  IMPORTED_LOCATION "${_zlib_lib}"
  INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${ZLIBNG_PREFIX}/include>"
)

add_dependencies(ZLIB::ZLIB zlibng_ep)

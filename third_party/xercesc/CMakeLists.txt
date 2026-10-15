include(ExternalProject)
include(GNUInstallDirs)

if (PACIFIC_USE_SUBMODULES)
  ExternalProject_Add(xercesc_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/xercesc/xercesc_submodule"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/xercesc-build"
    INSTALL_DIR "${XERECESC_PREFIX}"

    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
    BUILD_BYPRODUCTS
      "${XERECESC_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libxerces-c.so"
  )
else()
  ExternalProject_Add(xercesc_ep
    GIT_REPOSITORY https://github.com/apache/xerces-c.git
    GIT_TAG        xerces-3.3

    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/xercesc-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/xercesc-build"
    INSTALL_DIR "${XERECESC_PREFIX}"

    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
    BUILD_BYPRODUCTS
      "${XERECESC_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libxerces-c.so"
  )
endif()


file(MAKE_DIRECTORY "${XERECESC_PREFIX}/include")
file(MAKE_DIRECTORY "${XERECESC_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

add_library(XercesC::XercesC UNKNOWN IMPORTED GLOBAL)
set(_xercesc_xercesc "${XERECESC_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libxerces-c.so")
set_target_properties(XercesC::XercesC PROPERTIES
  IMPORTED_LOCATION "${_xercesc_xercesc}"
  INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${XERECESC_PREFIX}/include>"
)

add_dependencies(XercesC::XercesC xercesc_ep)
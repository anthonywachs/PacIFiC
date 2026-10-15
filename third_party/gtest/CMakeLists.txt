include(ExternalProject)
include(GNUInstallDirs)

if (PACIFIC_USE_SUBMODULES)
  ExternalProject_Add(gtest_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/gtest/gtest_submodule" 
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/gtest-build"
    INSTALL_DIR "${GTEST_PREFIX}"
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
else()
  ExternalProject_Add(gtest_ep
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.17.0

    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/gtest-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/gtest-build"
    INSTALL_DIR "${GTEST_PREFIX}"
  
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
endif()
file(MAKE_DIRECTORY "${GTEST_PREFIX}/include")
file(MAKE_DIRECTORY "${GTEST_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

add_library(GTest::gtest STATIC IMPORTED GLOBAL)
set_target_properties(GTest::gtest PROPERTIES
  IMPORTED_LOCATION "${GTEST_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${GTEST_PREFIX}/include"
)

add_dependencies(GTest::gtest gtest_ep)

add_library(GTest::gtest_main STATIC IMPORTED GLOBAL)
set_target_properties(GTest::gtest_main PROPERTIES
  IMPORTED_LOCATION "${GTEST_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_LINK_LIBRARIES "GTest::gtest"
  INTERFACE_INCLUDE_DIRECTORIES "${GTEST_PREFIX}/include"
)
add_dependencies(GTest::gtest_main gtest_ep)

list(PREPEND CMAKE_PREFIX_PATH "${GTEST_PREFIX}")

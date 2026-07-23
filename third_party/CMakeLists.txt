

option(BUILD_THIRD_PARTY_HDF5 "Build HDF5" OFF)
if (BUILD_THIRD_PARTY_HDF5)
    set(HDF5_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "HDF5 install prefix")
    add_subdirectory(hdf5)
endif()

option(BUILD_THIRD_PARTY_ZLIB "Build zlib" OFF)
if (BUILD_THIRD_PARTY_ZLIB)
    set(ZLIBNG_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for zlib-ng")
    add_subdirectory(zlib)
endif()

option(BUILD_THIRD_PARTY_XERCESC "Build XercesC" OFF)
if (BUILD_THIRD_PARTY_XERCESC)
    set(XERECESC_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Install prefix for XercesC")
    add_subdirectory(xercesc)
endif()

option(BUILD_THIRD_PARTY_PETSC "Build PETSc" OFF)
if(BUILD_THIRD_PARTY_PETSC)
    set(PETSC_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "PETSc install prefix")
    add_subdirectory(petsc)
endif()

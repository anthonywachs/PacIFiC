
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Build type: Debug, Release, RelWithDebInfo, or MinSizeRel." FORCE)
endif()

option(PACIFIC_BUILD_DOCS "Build documentation" OFF)
option(PACIFIC_BUILD_TESTS "Build unit & integration tests" ON)

option(PACIFIC_USE_MPI "Enable MPI-dependent targets" ON)
option(PACIFIC_USE_CUDA "Enable CUDA-dependent targets" OFF)
option(PACIFIC_USE_OPENMP "Enable OpenMP-dependent targets" ON)
option(PACIFIC_USE_GTEST "Enable GTest-depdendent targets" ON)

option(PACIFIC_BUILD_THIRD_PARTY_HDF5 "Build HDF5" OFF)
option(PACIFIC_BUILD_THIRD_PARTY_ZLIB "Build zlib" OFF)
option(PACIFIC_BUILD_THIRD_PARTY_XERCESC "Build XercesC" OFF)
option(PACIFIC_BUILD_THIRD_PARTY_GTEST "Build Google Test" OFF)
option(PACIFIC_BUILD_THIRD_PARTY_PETSC "Build PETSc" OFF)

option(PACIFIC_USE_SUBMODULES "Use git submodules for third-party dependencies where available" OFF)

option(PACIFIC_VERSION_USE_GIT "If no user override, compute version from git when available" ON)
option(PACIFIC_DEV_ENV "Create environment file" OFF) 


set(GRAINS3D_GPU_CUDA_STANDARD "17" CACHE STRING
  "Grains3DGPU CUDA C++ language standard")
set_property(CACHE GRAINS3D_GPU_CUDA_STANDARD PROPERTY STRINGS 14 17 20 23)
option(GRAINS3D_GPU_CUDA_STANDARD_REQUIRED
  "Grains3DGPU CUDA language standard required" ON)
option(GRAINS3D_GPU_CUDA_SEPARABLE_COMPILATION
  "Enable CUDA separable compilation" ON)
set(GRAINS3D_GPU_CUDA_ARCHITECTURES "75" CACHE STRING
  "CUDA architectures (semicolon-separated, e.g. 75;86)")
option(GRAINS3D_GPU_CUDA_RESOLVE_DEVICE_SYMBOLS
  "Enable CUDA device symbol resolution" ON)
set(GRAINS3D_GPU_CUDA_COMPILE_OPTIONS_LIST 
  "--extended-lambda;--expt-relaxed-constexpr" 
  CACHE STRING
  "CUDA compilation flags (semicolon-separated preferred; space-separated also accepted)"
)
if (GRAINS3D_GPU_CUDA_COMPILE_OPTIONS_LIST MATCHES " ")
  separate_arguments(GRAINS3D_GPU_CUDA_COMPILE_OPTIONS_LIST UNIX_COMMAND
    "${GRAINS3D_GPU_CUDA_COMPILE_OPTIONS_LIST}"
  )
endif()
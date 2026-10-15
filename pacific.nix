{
  lib,
  stdenv,
  fetchFromGitHub,
  cmake,
  pkg-config,
  gnumake,
  gcc,
  openmpi,
  hdf5-mpi,
  zlib,
  petsc,
  hypre,
  xercesc,
  makeWrapper,
  gtest,
  cudaToolkit ? null,
  cudaNvcc ? null,
}:

stdenv.mkDerivation {
  pname = "pacific";
  version = "0.0.1";

  # Build the current source tree (this repo checkout)
  src = ./.;

  nativeBuildInputs = [
    gcc
    cmake
    gnumake
    pkg-config
    makeWrapper
  ]
  ++ lib.optionals (cudaNvcc != null) [ cudaNvcc ];

  buildInputs = [
    gtest
    xercesc
    zlib
    (hdf5-mpi.override { mpi = openmpi; })
    openmpi
    petsc
  ]
  ++ lib.optionals (cudaToolkit != null) [ cudaToolkit ];

  propagatedBuildInputs = [
    gtest
    xercesc
    zlib
    (hdf5-mpi.override { mpi = openmpi; })
    openmpi
  ];

  cmakeFlags = [ 
    "-DCMAKE_BUILD_TYPE=Release"
    "-DUSE_SUBMODULES=ON"
    "-DOCTREE_BASILISK_PROVIDER=VENDORED"
  ]; 

  meta = with lib; {
    description = "PacIFiC multiphysics toolkit";
    homepage = "https://github.com/tortotubus/PacIFiC";
    platforms = platforms.linux;
    license = licenses.mit;
    mainProgram = "grains3d";
  };
}

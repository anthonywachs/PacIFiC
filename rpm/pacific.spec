
# PacIFiC Fedora Specfile (MPI flavors, CMake)
# ============================================

Name:           pacific
Version:        %{?version}
Release:        %{?release}
Summary:        PacIFiC toolkit (FSI + granular mechanics)
License:        MIT
URL:            https://gitlab.math.ubc.ca/pacific-devel-team/pacific
Source0:        pacific-%{version}.tar.gz

BuildRequires:  environment-modules
BuildRequires:  rpm-mpi-hooks
BuildRequires:  cmake-rpm-macros

BuildRequires:  openmpi-devel
BuildRequires:  hdf5-openmpi-devel
BuildRequires:  petsc-openmpi-devel

BuildRequires:  cmake >= 3.20
BuildRequires:  gcc-c++
BuildRequires:  make
BuildRequires:  xerces-c-devel
BuildRequires:  zlib-devel

BuildRequires: pacific-basilisk
 
%description
PacIFiC: high-performance fluid–structure interaction and granular mechanics toolkit.


%global __spec_build_shell /bin/bash

# ---------- Common (noarch) ----------
%package common
Summary: Common data files for PacIFiC
BuildArch: noarch

%description common

# ---------- OpenMPI flavor ----------
%package openmpi
Summary:        PacIFiC built with OpenMPI
Requires:       pacific-common = %{version}-%{release}
Requires:       openmpi
Requires:       hdf5-openmpi
Requires:       xerces-c
Requires:       zlib
Requires:       petsc-openmpi
Requires:       pacific-basilisk

%description openmpi
PacIFiC libraries and executables compiled using OpenMPI.

# ---------- OpenMPI flavor (devel) ----------
%package openmpi-devel
Summary:        Development files for PacIFiC (OpenMPI)
Requires:       %{name}-openmpi%{?_isa} = %{version}-%{release}
Requires:       openmpi-devel
Requires:       hdf5-openmpi-devel
Requires:       xerces-c-devel
Requires:       zlib-devel
Requires:       petsc-openmpi-devel

%description openmpi-devel
Headers and CMake package files for developing against PacIFiC built with OpenMPI.

%prep
%autosetup 

%build
%{_openmpi_load}
%global __cmake_builddir redhat-linux-build-openmpi
%if 0%{?rhel} < 10
export CC=mpicc CXX=mpicxx FC=mpifort
%endif 
%cmake \
  -DCMAKE_SKIP_RPATH:BOOL=ON \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF \
  -DCMAKE_INSTALL_RPATH:STRING= \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=OFF \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_MPI_INSTALL_LIBDIR:PATH=$MPI_LIB \
  -DCMAKE_MPI_INSTALL_BINDIR:PATH=$MPI_BIN \
  -DCMAKE_MPI_INSTALL_INCLUDEDIR:PATH=$MPI_INCLUDE \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
  -DPACIFIC_VERSION_SEMVER_OVERRIDE=%{?version} \
  -DPACIFIC_VERSION_FULL_OVERRIDE=%{?version}-%{?release} \
  -DOCTREE_BASILISK_PROVIDER=SYSTEM
%cmake_build 
%{_openmpi_unload}

%install
%{_openmpi_load} 
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component Grains3D_Runtime
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component Grains3D_Devel
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component Grains3D_Data
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component MAC_Runtime
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component MAC_Devel
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component FLUID_Runtime
DESTDIR=%{buildroot} cmake --install %{__cmake_builddir} --component Octree_Devel
%{_openmpi_unload}

%post openmpi -p /sbin/ldconfig
%postun openmpi -p /sbin/ldconfig

%files common
%license LICENSE
%doc README.md
%{_includedir}/Octree/*
%{_libdir}/cmake/Octree/*
%{_datadir}/Grains3D

%files openmpi
%{_libdir}/openmpi/bin/*
%{_libdir}/openmpi/lib/*.so*
%{_libdir}/openmpi/lib/cmake/Grains3D 
%{_libdir}/openmpi/lib/cmake/MAC

%files openmpi-devel
%{_includedir}/openmpi-x86_64/Grains3D
%{_includedir}/openmpi-x86_64/MAC

%changelog
* Thu Jan 15 2026 Conor Olive - %{version}-%{release}
- Fedora-style MPI flavor builds (openmpi/mpich)

#!/bin/bash
# Petsc installation script
# -------------------------

# Clean all
make allclean
if [ -d ${PETSC_ARCH} ]
then
  rm -rf ${PETSC_ARCH}
fi


# Create library links for Petsc
# ------------------------------
# Create libblas links for Petsc
LIBLAS_FOR_PETSC___=$(echo ${MACWORLD_BLAS_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
LIBLAS_FOR_PETSC___="-L${MACWORLD_BLAS_LIBDIR} $LIBLAS_FOR_PETSC___"
LIBATLAS_FOR_PETSC___=$(echo ${MACWORLD_ATLAS_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
LIBATLAS_FOR_PETSC___="-L${MACWORLD_ATLAS_LIBDIR} $LIBATLAS_FOR_PETSC___"
echo $'Blas/Atlas links for Petsc = ' ${LIBLAS_FOR_PETSC___} ${LIBATLAS_FOR_PETSC___}


# Create liblapack links for Petsc
LIBLAPACK_FOR_PETSC___=$(echo ${MACWORLD_LAPACK_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
LIBLAPACK_FOR_PETSC___="-L${MACWORLD_LAPACK_LIBDIR} $LIBLAPACK_FOR_PETSC___"
echo $'Lapack links for Petsc = ' ${LIBLAPACK_FOR_PETSC___}


# Create libmpi links for Petsc
LIBMPI_FOR_PETSC___=$(echo ${MACWORLD_MPI_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
LIBMPI_FOR_PETSC___="-L${MACWORLD_MPI_LIBDIR} $LIBMPI_FOR_PETSC___"
echo $'MPI links for Petsc = ' ${LIBMPI_FOR_PETSC___}


# Create libintel links for Petsc
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "Intel" ]]
then
  LIBINTEL_FOR_PETSC___=$(echo ${MACWORLD_INTEL_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
  LIBINTEL_FOR_PETSC___="-L${MACWORLD_INTEL_LIBDIR} $LIBINTEL_FOR_PETSC___"
  echo 'Intel links for Petsc =' $LIBINTEL_FOR_PETSC___
else
  LIBINTEL_FOR_PETSC___=""
  echo 'No Intel links for Petsc'
fi
export LIBINTEL_FOR_HYPRE___


# Create libgfortran links for HYPRE
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "GNU" ]]
then
  LIBGNU_FOR_PETSC___=$(echo ${MACWORLD_GFORTRAN_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
  LIBGNU_FOR_PETSC___="-L${MACWORLD_GFORTRAN_LIBDIR} $LIBGNU_FOR_PETSC___"
  echo 'GNU links for Petsc =' $LIBGNU_FOR_PETSC___
else
  LIBGNU_FOR_PETSC___=""
  echo 'No GNU links for Petsc'
fi
export LIBGNU_FOR_PETSC___

# ------------------------------


# Shared library ?
if [[ "${PETSC_SHARED_LIB}" == "yes" ]]
then
  boolshared="1"
else
  boolshared="0"
fi
echo 'Petsc as shared library =' ${boolshared}

# GNU compiler
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "GNU" ]]
then
  boolgnucompilers="1"
else
  boolgnucompilers="0"
fi
echo 'Petsc compiled with GNU =' ${boolgnucompilers}

# With HYPRE library as shared or not
if [[ "${HYPRE_SHARED_LIB}" == "yes" ]]
then
  hyprelibext="so"
else
  hyprelibext="a"
fi
echo 'HYPRE lib extension =' ${hyprelibext}


# Configure
# Options assumed that
# - HYPRE installed manually
# - if with MUMPS: MUMPS, Parmetis, PTScotch, Scalapack and Blacs downloaded from the Petsc website http://ftp.mcs.anl.gov/pub/petsc/externalpackages/
if [[ "${MACWORLD_PETSC_WITH_MUMPS}" == "1" ]]
then
  MUMPS_OPTS="--with-blacs=1 --download-blacs=yes --with-scalapack=1 --download-scalapack=yes --with-ptscotch=1 --download-ptscotch=yes --with-parmetis=1 --download-parmetis=yes --with-mumps=1 --download-mumps=yes"
else
  MUMPS_OPTS=""
fi
./configure --with-fortran --with-pic=1 --with-cmake=${PACIFIC_CMAKE} --CC="${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_C}" --CXX="${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_CXX}" --FC="${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_F90}" --COPTFLAGS="${PETSC_OPT_FLAGS}" --CXXOPTFLAGS="${PETSC_OPT_FLAGS}" --FOPTFLAGS="${PETSC_OPT_FLAGS}" --CFLAGS="-pthread" --CXXFLAGS="-pthread" --FFLAGS="-pthread" --with-debugging=0 --with-mpi=1 --with-shared-libraries=${boolshared} --with-mpi-compilers=1 --with-gnu-compilers=${boolgnucompilers} --with-mpiexec="${MACWORLD_MPI_BINDIR}/mpiexec" --with-blas-lapack-lib="${LIBLAS_FOR_PETSC___} ${LIBATLAS_FOR_PETSC___} ${LIBLAPACK_FOR_PETSC___}" --with-hypre=1  --with-hypre-include="${HYPRE_DIR}/${HYPRE_ARCH}/include" --with-hypre-lib="${HYPRE_DIR}/${HYPRE_ARCH}/lib/libHYPRE.${hyprelibext}" --with-mpi-include="[${MACWORLD_MPI_INCDIR},${MACWORLD_MPI_GFORTRAN_INCDIR}]" --with-mpi-lib="${LIBMPI_FOR_PETSC___}" --LDFLAGS="${LIBINTEL_FOR_PETSC___} ${LIBGNU_FOR_PETSC___}" ${MUMPS_OPTS}

# Compile and test
make
make test

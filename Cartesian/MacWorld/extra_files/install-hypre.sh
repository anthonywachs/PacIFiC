#!/bin/bash
# HYPRE installation script
#--------------------------

# Create library links for HYPRE
# ------------------------------
# Create libintel links for HYPRE
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "Intel" ]]
then
  LIBINTEL_FOR_HYPRE___=$(echo ${MACWORLD_INTEL_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
  LIBINTEL_FOR_HYPRE___="-L${MACWORLD_INTEL_LIBDIR} $LIBINTEL_FOR_HYPRE___"
  echo 'Intel links for Hypre =' $LIBINTEL_FOR_HYPRE___
else
  LIBINTEL_FOR_HYPRE___=""
  echo 'No Intel links for Hypre'
fi
export LIBINTEL_FOR_HYPRE___


# Create libgfortran links for HYPRE
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "GNU" ]]
then
  LIBGNU_FOR_HYPRE___=$(echo ${MACWORLD_GFORTRAN_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
  LIBGNU_FOR_HYPRE___="-L${MACWORLD_GFORTRAN_LIBDIR} $LIBGNU_FOR_HYPRE___"
  echo 'GNU links for Hypre =' $LIBGNU_FOR_HYPRE___
else
  LIBGNU_FOR_HYPRE___=""
  echo 'No GNU links for Hypre'
fi
export LIBGNU_FOR_HYPRE___
# ------------------------------


# Shared library ?
if [[ "${HYPRE_SHARED_LIB}" == "yes" ]]
then
  shared="--enable-shared"
fi


# Remove include & lib directory
if [ -d ${HYPRE_ARCH} ] 
then
  rm -rf ${HYPRE_ARCH}
fi


# Go to src directory
cd src

# Clean all
make clean
make distclean

# Configure
configure CC="${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_C}" CXX="${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_CXX}" F77="${MACWORLD_MPI_BINDIR}/${MACWORLD_MPI_F77}" CFLAGS="-fPIC ${HYPRE_OPT_FLAGS}" CXXFLAGS="-fPIC ${HYPRE_OPT_FLAGS}" FFLAGS="-fPIC ${HYPRE_OPT_FLAGS}" CPPFLAGS="-fPIC ${HYPRE_OPT_FLAGS}" --prefix=${HYPRE_DIR}/${HYPRE_ARCH} --libdir=${HYPRE_DIR}/${HYPRE_ARCH}/lib --includedir=${HYPRE_DIR}/${HYPRE_ARCH}/include --with-MPI-include=${MACWORLD_MPI_INCDIR} --with-MPI-libs="${MACWORLD_MPI_LIBS}" --with-MPI-lib-dirs=${MACWORLD_MPI_LIBDIR} --with-blas-lib-dirs="${MACWORLD_BLAS_LIBDIR} ${MACWORLD_ATLAS_LIBDIR}" --with-lapack-lib-dirs=${MACWORLD_LAPACK_LIBDIR} --with-blas-libs="${MACWORLD_BLAS_LIBS} ${MACWORLD_ATLAS_LIBS}" --with-lapack-libs="${MACWORLD_LAPACK_LIBS}" LDFLAGS="${LIBINTEL_FOR_HYPRE___} ${LIBGNU_FOR_HYPRE___}" $shared

# Compile and install
make install | tee -a compil.log
mv compil.log ../${HYPRE_ARCH}/
cp config.log ../${HYPRE_ARCH}/

# Go back to root directory
cd ../

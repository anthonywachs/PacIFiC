# MacWorld root
export MACWORLD_ROOT=${PACIFIC_HOME}/Cartesian/MacWorld
export MACWORLD_BITS_DEFAULT=${PACIFIC_BITS_DEFAULT}
export MACWORLD_BITS_EXT=${PACIFIC_BITS_EXT}

# MPI
export MACWORLD_MPI_ROOT=${PACIFIC_MPI_ROOT}
export MACWORLD_MPI_DISTRIB=${PACIFIC_MPI_DISTRIB}
export MACWORLD_MPI_VERSION=${PACIFIC_MPI_VERSION}
export MACWORLD_MPI_INCDIR=${PACIFIC_MPI_INCDIR}
export MACWORLD_MPI_GFORTRAN_INCDIR=${PACIFIC_MPI_GFORTRAN_INCDIR}
export MACWORLD_MPI_BINDIR=${PACIFIC_MPI_BINDIR}
export MACWORLD_MPI_LIBDIR=${PACIFIC_MPI_LIBDIR}
export MACWORLD_MPI_C=${PACIFIC_MPI_C}
export MACWORLD_MPI_CXX=${PACIFIC_MPI_CXX}
export MACWORLD_MPI_F77=${PACIFIC_MPI_F77}
export MACWORLD_MPI_F90=${PACIFIC_MPI_F90}
export MACWORLD_MPI_LIBS=${PACIFIC_MPI_LIBS}
export MACWORLD_MPI_CPPLIBS=${PACIFIC_MPI_CPPLIBS}
export MACWORLD_MPI_CFLIBS=${PACIFIC_MPI_CFLIBS}

# Serial compiler
export MACWORLD_SERCOMPIL_ENV=${PACIFIC_SERCOMPIL_ENV}
export MACWORLD_SERCOMPIL_VERSION=${PACIFIC_SERCOMPIL_VERSION}
export MACWORLD_OPT_FLAGS=${PACIFIC_OPT_FLAGS}

# MacWorld full extension
if [ "${MACWORLD_BITS_DEFAULT}" == "64" ]
then
  MACWORLD_PRE="${MACWORLD_BITS_DEFAULT}-"
fi
export MACWORLD_FULL_EXT=${MACWORLD_PRE}${MACWORLD_MPI_DISTRIB}-${MACWORLD_MPI_VERSION}-${MACWORLD_SERCOMPIL_ENV}-${MACWORLD_SERCOMPIL_VERSION}

# Low level libraries for MacWorld for PeliPACK
export MACWORLD_BLAS_LIBDIR=${PACIFIC_BLAS_LIBDIR}
export MACWORLD_BLAS_LIBS=${PACIFIC_BLAS_LIBS}
export MACWORLD_ATLAS_LIBDIR=${PACIFIC_ATLAS_LIBDIR}
export MACWORLD_ATLAS_LIBS=${PACIFIC_ATLAS_LIBS}
export MACWORLD_LAPACK_LIBDIR=${PACIFIC_LAPACK_LIBDIR}
export MACWORLD_LAPACK_LIBS=${PACIFIC_LAPACK_LIBS}
export MACWORLD_GFORTRAN_LIBDIR=${PACIFIC_GFORTRAN_LIBDIR}
export MACWORLD_GFORTRAN_LIBS=${PACIFIC_GFORTRAN_LIBS}
export MACWORLD_INTEL_LIBDIR=${PACIFIC_INTEL_LIBDIR}
export MACWORLD_INTEL_LIBS=${PACIFIC_INTEL_LIBS}
export MACWORLD_M_LIBDIR=${PACIFIC_M_LIBDIR}
export MACWORLD_Z_DIR=${PACIFIC_Z_DIR}
export MACWORLD_Z_INCDIR=${PACIFIC_Z_INCDIR}
export MACWORLD_Z_LIBDIR=${PACIFIC_Z_LIBDIR}
export MACWORLD_X11_DIR=${PACIFIC_X11_DIR}
export MACWORLD_X11_INCDIR=${PACIFIC_X11_INCDIR}
export MACWORLD_X11_LIBDIR=${PACIFIC_X11_LIBDIR}

# Other
export MACWORLD_PETSC_WITH_MUMPS=${PACIFIC_PETSC_WITH_MUMPS}

echo -e '\033[93m*** MacWorld shell variables\033[0m'
echo -e '\033[93mMACWORLD_ROOT\033[0m =' $MACWORLD_ROOT
echo -e '\033[93mMACWORLD_BITS_DEFAULT\033[0m =' $MACWORLD_BITS_DEFAULT
echo -e '\033[93mMACWORLD_BITS_EXT\033[0m =' $MACWORLD_BITS_EXT
echo -e '\033[93mMACWORLD_MPI_ROOT\033[0m =' $MACWORLD_MPI_ROOT
echo -e '\033[93mMACWORLD_MPI_DISTRIB\033[0m =' $MACWORLD_MPI_DISTRIB
echo -e '\033[93mMACWORLD_MPI_VERSION\033[0m =' $MACWORLD_MPI_VERSION
echo -e '\033[93mMACWORLD_MPI_INCDIR\033[0m =' $MACWORLD_MPI_INCDIR
echo -e '\033[93mMACWORLD_MPI_GFORTRAN_INCDIR\033[0m =' $MACWORLD_MPI_GFORTRAN_INCDIR
echo -e '\033[93mMACWORLD_MPI_BINDIR\033[0m =' $MACWORLD_MPI_BINDIR
echo -e '\033[93mMACWORLD_MPI_LIBDIR\033[0m =' $MACWORLD_MPI_LIBDIR
echo -e '\033[93mMACWORLD_FULL_EXT\033[0m =' $MACWORLD_FULL_EXT
echo -e '\033[93mMACWORLD_MPI_C\033[0m =' $MACWORLD_MPI_C
echo -e '\033[93mMACWORLD_MPI_CXX\033[0m =' $MACWORLD_MPI_CXX
echo -e '\033[93mMACWORLD_MPI_F77\033[0m =' $MACWORLD_MPI_F77
echo -e '\033[93mMACWORLD_MPI_F90\033[0m =' $MACWORLD_MPI_F90
echo -e '\033[93mMACWORLD_MPI_LIBS\033[0m =' $MACWORLD_MPI_LIBS
echo -e '\033[93mMACWORLD_MPI_CPPLIBS\033[0m =' $MACWORLD_MPI_CPPLIBS
echo -e '\033[93mMACWORLD_MPI_CFLIBS\033[0m =' $MACWORLD_MPI_CFLIBS
echo -e '\033[93mMACWORLD_BLAS_LIBDIR\033[0m =' $MACWORLD_BLAS_LIBDIR
echo -e '\033[93mMACWORLD_BLAS_LIBS\033[0m =' $MACWORLD_BLAS_LIBS
echo -e '\033[93mMACWORLD_ATLAS_LIBDIR\033[0m =' $MACWORLD_ATLAS_LIBDIR
echo -e '\033[93mMACWORLD_ATLAS_LIBS\033[0m =' $MACWORLD_ATLAS_LIBS
echo -e '\033[93mMACWORLD_LAPACK_LIBDIR\033[0m =' $MACWORLD_LAPACK_LIBDIR
echo -e '\033[93mMACWORLD_LAPACK_LIBS\033[0m =' $MACWORLD_LAPACK_LIBS
echo -e '\033[93mMACWORLD_GFORTRAN_LIBDIR\033[0m =' $MACWORLD_GFORTRAN_LIBDIR
echo -e '\033[93mMACWORLD_GFORTRAN_LIBS\033[0m =' $MACWORLD_GFORTRAN_LIBS
echo -e '\033[93mMACWORLD_INTEL_LIBDIR\033[0m =' $MACWORLD_INTEL_LIBDIR
echo -e '\033[93mMACWORLD_INTEL_LIBS\033[0m =' $MACWORLD_INTEL_LIBS
echo -e '\033[93mMACWORLD_M_LIBDIR\033[0m =' $MACWORLD_M_LIBDIR
echo -e '\033[93mMACWORLD_Z_DIR\033[0m =' $MACWORLD_Z_DIR
echo -e '\033[93mMACWORLD_Z_INCDIR\033[0m =' $MACWORLD_Z_INCDIR
echo -e '\033[93mMACWORLD_Z_LIBDIR\033[0m =' $MACWORLD_Z_LIBDIR
echo -e '\033[93mMACWORLD_X11_DIR\033[0m =' $MACWORLD_X11_DIR
echo -e '\033[93mMACWORLD_X11_INCDIR\033[0m =' $MACWORLD_X11_INCDIR
echo -e '\033[93mMACWORLD_X11_LIBDIR\033[0m =' $MACWORLD_X11_LIBDIR
echo -e '\033[93mMACWORLD_SERCOMPIL_ENV\033[0m =' $MACWORLD_SERCOMPIL_ENV
echo -e '\033[93mMACWORLD_SERCOMPIL_VERSION\033[0m =' $MACWORLD_SERCOMPIL_VERSION
echo -e '\033[93mMACWORLD_OPT_FLAGS\033[0m =' $MACWORLD_OPT_FLAGS
echo -e '\033[93mMACWORLD_PETSC_WITH_MUMPS\033[0m =' $MACWORLD_PETSC_WITH_MUMPS
echo -e '  '


# HYPRE / Petsc
echo -e '\033[96m*** Petsc/Hypre shell variables\033[0m'
source ${MACWORLD_ROOT}/hypre-2.24.0/hypre.env.sh
source ${MACWORLD_ROOT}/petsc-3.17.5/petsc.env.sh
echo -e '  '


# MAC
source ${MACWORLD_ROOT}/MAC/Env/init.sh
if [[ ${PACIFIC_AUTO_CONFIG} -eq 1 ]]
then
  echo -e '\033[33mUsing MAC template Linux and extra-Linux makefiles\033[0m'
  cp ${MACWORLD_ROOT}/MAC/etc/Linux_template.mak ${MACWORLD_ROOT}/MAC/etc/${MAC_ARCH}.mak
  cp ${MACWORLD_ROOT}/MAC/etc/extra-Linux_template.mak ${MACWORLD_ROOT}/MAC/etc/extra-${MAC_ARCH}.mak
fi

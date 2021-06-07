#
# PacIFiC 64 bits + OpenMPI + GNU
#

# Modules
module purge
module load StdEnv/2018.3
module load gcc/7.3.0
module load openmpi/2.1.1
module load imkl/11.3.4.258


# PacIFiC home
export PACIFIC_HOME=/home/awachs/PacIFiC-Git/pacific
export PACIFIC_EXE_SCRIPTS=${PACIFIC_HOME}/ExeScripts
export PATH="${PACIFIC_EXE_SCRIPTS}:${PATH}"
export PACIFIC_BITS_DEFAULT="64"
export PACIFIC_BITS_EXT="64"
echo -e '\033[94m*** PacIFiC shell variables\033[0m'
echo -e '\033[94mPACIFIC_HOME\033[0m =' $PACIFIC_HOME
echo -e '\033[94mPACIFIC_BITS_DEFAULT\033[0m =' $PACIFIC_BITS_DEFAULT
echo -e '\033[94mPACIFIC_BITS_EXT\033[0m =' $PACIFIC_BITS_EXT
echo -e '  '


# MPI
export PACIFIC_MPI_ROOT=/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/openmpi/2.1.1
export PACIFIC_MPI_DISTRIB=OpenMPI
export PACIFIC_MPI_VERSION="2.1.1"
export PACIFIC_MPI_INCDIR="${PACIFIC_MPI_ROOT}/include"
export PACIFIC_MPI_GFORTRAN_INCDIR="${PACIFIC_MPI_ROOT}/include"
export PACIFIC_MPI_BINDIR="${PACIFIC_MPI_ROOT}/bin"
export PATH="${PACIFIC_MPI_BINDIR}:${PATH}"
export PACIFIC_MPI_LIBDIR="${PACIFIC_MPI_ROOT}/lib"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PACIFIC_MPI_LIBDIR}"
export PACIFIC_MPI_C="mpicc"
export PACIFIC_MPI_CXX="mpic++"
export PACIFIC_MPI_F77="mpifort"
export PACIFIC_MPI_F90="mpifort" 
export PACIFIC_MPI_LIBS="mpi mpi_cxx mpi_mpifh"
export PACIFIC_MPI_CPPLIBS="mpi mpi_cxx"
export PACIFIC_MPI_CFLIBS="mpi mpi_mpifh"
echo -e '\033[32m*** MPI shell variables\033[0m'
echo -e '\033[32mPACIFIC_MPI_ROOT\033[0m =' $PACIFIC_MPI_ROOT
echo -e '\033[32mPACIFIC_MPI_DISTRIB\033[0m =' $PACIFIC_MPI_DISTRIB
echo -e '\033[32mPACIFIC_MPI_VERSION\033[0m =' $PACIFIC_MPI_VERSION
echo -e '\033[32mPACIFIC_MPI_INCDIR\033[0m =' $PACIFIC_MPI_INCDIR
echo -e '\033[32mPACIFIC_MPI_GFORTRAN_INCDIR\033[0m =' $PACIFIC_MPI_GFORTRAN_INCDIR
echo -e '\033[32mPACIFIC_MPI_BINDIR\033[0m =' $PACIFIC_MPI_BINDIR
echo -e '\033[32mPACIFIC_MPI_LIBDIR\033[0m =' $PACIFIC_MPI_LIBDIR
echo -e '\033[32mPACIFIC_MPI_C\033[0m =' $PACIFIC_MPI_C
echo -e '\033[32mPACIFIC_MPI_CXX\033[0m =' $PACIFIC_MPI_CXX
echo -e '\033[32mPACIFIC_MPI_F77\033[0m =' $PACIFIC_MPI_F77
echo -e '\033[32mPACIFIC_MPI_F90\033[0m =' $PACIFIC_MPI_F90
echo -e '\033[32mPACIFIC_MPI_LIBS\033[0m =' $PACIFIC_MPI_LIBS
echo -e '\033[32mPACIFIC_MPI_CPPLIBS\033[0m =' $PACIFIC_MPI_CPPLIBS
echo -e '\033[32mPACIFIC_MPI_CFLIBS\033[0m =' $PACIFIC_MPI_CFLIBS
echo -e '  '


# Serial compiler and low level librairies
export PACIFIC_SERCOMPIL_ENV="GNU"
export PACIFIC_SERCOMPIL_VERSION="7.3.0"
if [[ "${PACIFIC_SERCOMPIL_ENV}" == "GNU" ]] 
then
  PACIFIC_SERCOMPIL_C=$(which gcc)
  PACIFIC_SERCOMPIL_CPP=$(which g++)  
else if [[ "${PACIFIC_SERCOMPIL_ENV}" == "Intel" ]]  
  then
    PACIFIC_SERCOMPIL_C=$(which icc)
    PACIFIC_SERCOMPIL_CPP=$(which icpc)    
  else
    PACIFIC_SERCOMPIL_C="C compiler undefined"
    PACIFIC_SERCOMPIL_CPP="C++ compiler undefined"    
  fi
fi
export PACIFIC_SERCOMPIL_C
export PACIFIC_SERCOMPIL_CPP
export PACIFIC_OPT_FLAGS="-O3"
export PACIFIC_BLAS_LIBDIR=/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/11.3.4.258/mkl/lib/intel64
export PACIFIC_BLAS_LIBS="mkl_blas95_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_gf_lp64"
export PACIFIC_ATLAS_LIBDIR=/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/11.3.4.258/mkl/lib/intel64
export PACIFIC_ATLAS_LIBS="mkl_blas95_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_gf_lp64"
export PACIFIC_LAPACK_LIBDIR=/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/11.3.4.258/mkl/lib/intel64
export PACIFIC_LAPACK_LIBS="mkl_lapack95_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_gf_lp64"
export PACIFIC_GFORTRAN_LIBDIR=/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/gcc-7.3.0/lib64
export PACIFIC_GFORTRAN_LIBS="gfortran"
export PACIFIC_INTEL_LIBDIR=""
export PACIFIC_INTEL_LIBS=""
export PACIFIC_M_LIBDIR=/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib64
export PACIFIC_Z_DIR=/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09
export PACIFIC_Z_INCDIR="${PACIFIC_Z_DIR}/include"
export PACIFIC_Z_LIBDIR="${PACIFIC_Z_DIR}/lib"
export PACIFIC_X11_DIR=/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09
export PACIFIC_X11_INCDIR="${PACIFIC_X11_DIR}/include/X11"
export PACIFIC_X11_LIBDIR="${PACIFIC_X11_DIR}/lib"

echo -e '\033[32m*** Serial compiler and low level librairies shell variables\033[0m'
echo -e '\033[32mPACIFIC_SERCOMPIL_ENV\033[0m =' $PACIFIC_SERCOMPIL_ENV
echo -e '\033[32mPACIFIC_SERCOMPIL_VERSION\033[0m =' $PACIFIC_SERCOMPIL_VERSION
echo -e '\033[32mPACIFIC_SERCOMPIL_C\033[0m =' $PACIFIC_SERCOMPIL_C
echo -e '\033[32mPACIFIC_SERCOMPIL_CPP\033[0m =' $PACIFIC_SERCOMPIL_CPP
echo -e '\033[32mPACIFIC_OPT_FLAGS\033[0m =' $PACIFIC_OPT_FLAGS
echo -e '\033[32mPACIFIC_BLAS_LIBDIR\033[0m =' $PACIFIC_BLAS_LIBDIR
echo -e '\033[32mPACIFIC_BLAS_LIBS\033[0m =' $PACIFIC_BLAS_LIBS
echo -e '\033[32mPACIFIC_ATLAS_LIBDIR\033[0m =' $PACIFIC_ATLAS_LIBDIR
echo -e '\033[32mPACIFIC_ATLAS_LIBS\033[0m =' $PACIFIC_ATLAS_LIBS
echo -e '\033[32mPACIFIC_LAPACK_LIBDIR\033[0m =' $PACIFIC_LAPACK_LIBDIR
echo -e '\033[32mPACIFIC_LAPACK_LIBS\033[0m =' $PACIFIC_LAPACK_LIBS
echo -e '\033[32mPACIFIC_GFORTRAN_LIBDIR\033[0m =' $PACIFIC_GFORTRAN_LIBDIR
echo -e '\033[32mPACIFIC_GFORTRAN_LIBS\033[0m =' $PACIFIC_GFORTRAN_LIBS
echo -e '\033[32mPACIFIC_INTEL_LIBDIR\033[0m =' $PACIFIC_INTEL_LIBDIR
echo -e '\033[32mPACIFIC_INTEL_LIBS\033[0m =' $PACIFIC_INTEL_LIBS
echo -e '\033[32mPACIFIC_M_LIBDIR\033[0m =' $PACIFIC_M_LIBDIR
echo -e '\033[32mPACIFIC_Z_DIR\033[0m =' $PACIFIC_Z_DIR
echo -e '\033[32mPACIFIC_Z_INCDIR\033[0m =' $PACIFIC_Z_INCDIR
echo -e '\033[32mPACIFIC_Z_LIBDIR\033[0m =' $PACIFIC_Z_LIBDIR
echo -e '\033[32mPACIFIC_X11_DIR\033[0m =' $PACIFIC_X11_DIR
echo -e '\033[32mPACIFIC_X11_INCDIR\033[0m =' $PACIFIC_X11_INCDIR
echo -e '\033[32mPACIFIC_X11_LIBDIR\033[0m =' $PACIFIC_X11_LIBDIR
echo -e '  '


# Other
export PACIFIC_AUTO_CONFIG=1
export PACIFIC_PETSC_WITH_MUMPS=0
export PACIFIC_CMAKE=$(which cmake)
echo -e '\033[37m*** Other variables\033[0m'
echo -e '\033[37mPACIFIC_AUTO_CONFIG\033[0m =' $PACIFIC_AUTO_CONFIG
echo -e '\033[37mPACIFIC_PETSC_WITH_MUMPS\033[0m =' $PACIFIC_PETSC_WITH_MUMPS
echo -e '\033[37mPACIFIC_CMAKE\033[0m =' $PACIFIC_CMAKE
echo -e '  '


# Grains
if [[ ${PACIFIC_AUTO_CONFIG} -eq 1 ]]
then
  echo -e '\033[31mUsing grains_env_template.env.sh as env file\033[0m'
  cp ${PACIFIC_HOME}/GRAINS/Env/grains_env_template.env.sh ${PACIFIC_HOME}/GRAINS/Env/grains-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
fi
source ${PACIFIC_HOME}/GRAINS/Env/grains-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
echo -e '  '


# MacWorld
if [[ ${PACIFIC_AUTO_CONFIG} -eq 1 ]]
then
  echo -e '\033[93mUsing macworld_env_template.env.sh as env file\033[0m'
  cp ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld_env_template.env.sh ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
fi
source ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
echo -e '  '


# Fluid MAC
if [[ ${PACIFIC_AUTO_CONFIG} -eq 1 ]]
then
  echo -e '\033[90mUsing fluid_env_template.env.sh as env file\033[0m'
  cp ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid_env_template.env.sh ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
fi
source ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
echo -e '  '


# Basilisk
source ${PACIFIC_HOME}/Octree/Env/octree.env.sh

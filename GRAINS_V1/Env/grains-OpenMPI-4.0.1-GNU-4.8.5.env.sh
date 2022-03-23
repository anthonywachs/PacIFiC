# Definition
# Grains
export GRAINS_HOME=${HOME}/Multiphase/PacNewGrains-Git/pacific/GRAINS_V1
export GRAINS_BITS_DEFAULT=64
export GRAINS_SERCOMPIL_ENV="GNU"
export GRAINS_SERCOMPIL_VERSION="4.8.5"
export GRAINS_MPICCC=mpicxx
export GRAINS_CPP=g++
export GRAINS_C=gcc
export GRAINS_COMPIL_OPT="-O3"

# MPI
export GRAINS_MPI_ROOT=/home/wachs/local/openmpi-4.0.1
export GRAINS_MPI_DISTRIB=OpenMPI
export GRAINS_MPI_VERSION="4.0.1"
export GRAINS_MPI_INCDIR="${GRAINS_MPI_ROOT}/include"
export GRAINS_MPI_BINDIR="${GRAINS_MPI_ROOT}/bin"
export GRAINS_MPI_LIBDIR="${GRAINS_MPI_ROOT}/lib"
export GRAINS_MPI_LIBS="mpi"
# End definition


# Other environment variables
source $GRAINS_HOME/Env/grains_default.env.sh

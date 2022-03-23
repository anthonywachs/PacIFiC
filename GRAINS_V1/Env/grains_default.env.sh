# Bits extension
GRAINS_BITS_EXT=""
if [ $GRAINS_BITS_DEFAULT = "64" ] 
then
  GRAINS_BITS_EXT=64
fi
export GRAINS_BITS_EXT
# End Bits extension

# MPI linking libs
export LIBMPI_FOR_GRAINS___=$(echo -e ${GRAINS_MPI_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
# End MPI linking libs

# Xerces definition
export XERCESCROOT=${GRAINS_HOME}/XERCES-2.8.0
export XERCESC_ROOT=${XERCESCROOT}
# End Xerces definition

# Other definitions
export GRAINS_ROOT=${GRAINS_HOME}/Grains
# End other definitions

# Full extension
export GRAINS_FULL_EXT=${GRAINS_BITS_EXT}-${GRAINS_MPI_DISTRIB}-${GRAINS_MPI_VERSION}-${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION}
# End Full extension

# LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_ROOT}/lib${GRAINS_FULL_EXT}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_MPI_LIBDIR}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XERCESC_ROOT}/lib${GRAINS_BITS_EXT}-${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION}
# End LD_LIBRARY_PATH

# Display
echo -e '\033[31mGRAINS_HOME\033[0m =' $GRAINS_HOME
echo -e '\033[31mGRAINS_ROOT\033[0m =' $GRAINS_ROOT
echo -e '\033[31mGRAINS_BITS_DEFAULT\033[0m =' $GRAINS_BITS_DEFAULT
echo -e '\033[31mGRAINS_BITS_EXT\033[0m =' $GRAINS_BITS_EXT
echo -e '\033[31mGRAINS_SERCOMPIL_ENV\033[0m =' $GRAINS_SERCOMPIL_ENV
echo -e '\033[31mGRAINS_SERCOMPIL_VERSION\033[0m =' $GRAINS_SERCOMPIL_VERSION
echo -e '\033[31mGRAINS_CPP\033[0m =' $GRAINS_CPP
echo -e '\033[31mGRAINS_C\033[0m =' $GRAINS_C
echo -e '\033[31mGRAINS_COMPIL_OPT\033[0m =' $GRAINS_COMPIL_OPT
echo -e '\033[31mGRAINS_MPI_ROOT\033[0m =' $GRAINS_MPI_ROOT
echo -e '\033[31mGRAINS_MPI_DISTRIB\033[0m =' $GRAINS_MPI_DISTRIB
echo -e '\033[31mGRAINS_MPI_VERSION\033[0m =' $GRAINS_MPI_VERSION
echo -e '\033[31mGRAINS_MPI_INCDIR\033[0m =' $GRAINS_MPI_INCDIR
echo -e '\033[31mGRAINS_MPI_BINDIR\033[0m =' $GRAINS_MPI_BINDIR
echo -e '\033[31mGRAINS_MPI_LIBDIR\033[0m =' $GRAINS_MPI_LIBDIR
echo -e '\033[31mGRAINS_MPI_LIBS\033[0m =' $GRAINS_MPI_LIBS
echo -e '\033[31mLIBMPI_FOR_GRAINS___\033[0m =' $LIBMPI_FOR_GRAINS___
echo -e '\033[31mGRAINS_MPICCC\033[0m =' $GRAINS_MPICCC
echo -e '\033[31mGRAINS_FULL_EXT\033[0m =' $GRAINS_FULL_EXT
echo -e '\033[31mXERCESCROOT & XERCESC_ROOT\033[0m =' $XERCESCROOT '&' $XERCESC_ROOT


# General compilation flags - not displayed
if [ $GRAINS_CPP = "g++" ] 
then
  GCCFLAGS="-pedantic -W -Wno-long-long -Wno-ctor-dtor-privacy -Wno-unused-parameter -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11 "
fi
export GRAINS_GENCCFLAGS="-m${GRAINS_BITS_DEFAULT} ${GRAINS_COMPIL_OPT} -fPIC -Wall -Wwrite-strings -Wconversion -Wshadow -Wno-deprecated -Wno-comment ${GCCFLAGS}"
export GRAINS_MPICCFLAGS="-DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX"


# System include directories to generate dependencies 
GRAINS_SYSINC="${GRAINS_MPI_INCDIR} ${XERCESCROOT}/include"
# Next line depends on the operating system. Not mandatory to change, 
# unless you do not want makedepend to return 1000s of harmless warnings 
GRAINS_SYSINC="${GRAINS_SYSINC} /usr/include/c++/4.8.5  /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include /usr/include/c++/4.8.2/x86_64-redhat-linux"
export GRAINS_SYSINC

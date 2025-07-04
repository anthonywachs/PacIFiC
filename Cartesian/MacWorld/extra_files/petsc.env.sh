export PETSC_VERSION=3.2.0
export PETSC_VERSION_PATCH=3.2.0-p7
export PETSC_DIR=${MACWORLD_ROOT}/petsc-${PETSC_VERSION_PATCH}
export PETSC_HOME=${PETSC_DIR}
export PETSC_OPT_FLAGS=${MACWORLD_OPT_FLAGS}
export PETSC_SHARED_LIB="yes"
if [[ "${MACWORLD_PETSC_WITH_MUMPS}" == "1" ]]
then
  export MUMPS_VERSION=4.10.0
  export PETSC_ARCH=Linux-${MACWORLD_FULL_EXT}-HYPRE-${HYPRE_VERSION}-MUMPS-${MUMPS_VERSION}
else
  export MUMPS_VERSION=""
  export PETSC_ARCH=Linux-${MACWORLD_FULL_EXT}-HYPRE-${HYPRE_VERSION}
fi
export MUMPS_VERSION
export PETSC_ARCH

# Display
if [[ "${MACWORLD_PETSC_WITH_MUMPS}" == "1" ]]
then
  echo -e '\033[96mMUMPS_VERSION\033[0m =' $MUMPS_VERSION
fi
echo -e '\033[96mPETSC_VERSION\033[0m =' $PETSC_VERSION
echo -e '\033[96mPETSC_VERSION_PATCH\033[0m =' $PETSC_VERSION_PATCH
echo -e '\033[96mPETSC_DIR\033[0m =' $PETSC_DIR
echo -e '\033[96mPETSC_HOME\033[0m =' $PETSC_HOME
echo -e '\033[96mPETSC_OPT_FLAGS\033[0m =' $PETSC_OPT_FLAGS
echo -e '\033[96mPETSC_SHARED_LIB\033[0m =' $PETSC_SHARED_LIB
echo -e '\033[96mPETSC_ARCH\033[0m =' $PETSC_ARCH

# LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PETSC_DIR}/${PETSC_ARCH}/lib"

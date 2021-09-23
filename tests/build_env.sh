cd ../Env
sed -i "s|TO_BE_OVERWRITTEN|$(echo $PWD)|" ./CI-RUNNER-PacIFiC-OpenMPI-2.1.1-GNU-8.2.1.env.sh
source CI-RUNNER-PacIFiC-OpenMPI-2.1.1-GNU-8.2.1.env.sh
cp ${PACIFIC_HOME}/GRAINS/Env/grains_env_template.env.sh ${PACIFIC_HOME}/GRAINS/Env/grains-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
cp ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld_env_template.env.sh ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
cp ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid_env_template.env.sh ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
cd ..

# Modification of the environment files for the specific architecture of r8k1-wachs1.math.ubc.ca
sed -i 's|${GRAINS_HOME}/XERCES-2.8.0|/home/gitlab-runner/dependencies/XERCES-2.8.0|' GRAINS/Env/grains_default.env.sh
sed -i 's|${GRAINS_HOME}/XERCES-2.8.0|/home/gitlab-runner/dependencies/XERCES-2.8.0|' GRAINS/Env/grains_default.env.csh
sed -i 's|source ${MACWORLD_ROOT}/hypre-2.10.1/hypre.env.sh|source ${MACWORLD_ROOT}/extra_files/hypre.env.sh|' ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
sed -i 's|source ${MACWORLD_ROOT}/petsc-3.2.0-p7/petsc.env.sh|source ${MACWORLD_ROOT}/extra_files/petsc.env.sh|' ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
sed -i 's|HYPRE_DIR=${MACWORLD_ROOT}/hypre-2.10.1|HYPRE_DIR=/home/gitlab-runner/dependencies/hypre-2.10.1|' ${PACIFIC_HOME}/Cartesian/MacWorld/extra_files/hypre.env.sh
sed -i 's|PETSC_DIR=${MACWORLD_ROOT}/petsc-${PETSC_VERSION_PATCH}|PETSC_DIR=/home/gitlab-runner/dependencies/petsc-3.2.0-p7|' ${PACIFIC_HOME}/Cartesian/MacWorld/extra_files/petsc.env.sh
sed -i 's|export BASILISK=${OCTREE_HOME}/basilisk/src|export BASILISK=/home/gitlab-runner/dependencies/basilisk/src|' ${PACIFIC_HOME}/Octree/Env/octree.env.sh
sed -i 's|export BASILISK_MESA_GLU_INSTALL=|export BASILISK_MESA_GLU_INSTALL=/home/gitlab-runner/dependencies/local|' ${PACIFIC_HOME}/Octree/Env/octree.env.sh

mv Env/CI-RUNNER-PacIFiC-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh Env/PacIFiC-CI-RUNNER-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
source Env/PacIFiC-CI-RUNNER-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh

mv ${MAC_HOME}/etc/Linux_template.mak ${MAC_HOME}/etc/Linux-${MAC_FULL_EXT}.mak
mv ${MAC_HOME}/etc/extra-Linux_template.mak ${MAC_HOME}/etc/extra-Linux-${MAC_FULL_EXT}.mak

if [ ! -L ${BASILISK}/eulerian_caps ]
then
  ln -s ${OCTREE_HOME}/eulerian_caps ${BASILISK}/eulerian_caps
fi

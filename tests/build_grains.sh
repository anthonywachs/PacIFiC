source ../Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh

# Compilation of Grains3D
cd $GRAINS_HOME
./makeARCH create ; make update ; make dtd

# Check if the compilation was successful
test -f ${GRAINS_HOME}/Main/bin${GRAINS_FULL_EXT}/grains ; echo $? >> ${PACIFIC_HOME}/tests/build_success.txt

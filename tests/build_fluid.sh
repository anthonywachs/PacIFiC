source ../Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh

# Compilation of FLUID
cd ${PACIFIC_HOME}/Cartesian/FLUID
./compil
test -f ${PACIFIC_HOME}/Cartesian/FLUID/lib/Linux-${MAC_FULL_EXT}/exe0 ; echo $? >> ${PACIFIC_HOME}/tests/build_success.txt
test -f ${PACIFIC_HOME}/Cartesian/FLUID/lib/Linux-${MAC_FULL_EXT}/exe2 ; echo $? >> ${PACIFIC_HOME}/tests/build_success.txt

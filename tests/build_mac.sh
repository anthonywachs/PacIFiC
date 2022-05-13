source ../Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh

# Compilation of MacWorld
cd $MACWORLD_ROOT/MAC
./install-mac.sh

# Check if compilation was successful
test -f ${MAC_HOME}/lib/Linux-${MAC_FULL_EXT}/libmac0.so ; echo $? >> ${PACIFIC_HOME}/tests/build_success.txt
test -f ${MAC_HOME}/lib/Linux-${MAC_FULL_EXT}/libmac2.so ; echo $? >> ${PACIFIC_HOME}/tests/build_success.txt

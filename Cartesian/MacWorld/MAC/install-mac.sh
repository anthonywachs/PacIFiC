#!/bin/bash
# MAC installation script
# ----------------------------

## Keep previous architecture-based installations
__BACKUP_DIR="___BackUp_MACLib"
if [ -d lib ]
  then
    mkdir ${__BACKUP_DIR}
    cd lib/
    for dir in *
      do
	if [ "$dir" != "Linux-${MAC_FULL_EXT}" ]
          then
            cp -r $dir ../${__BACKUP_DIR}/
        fi
      done
    cd ../  
fi  


# Clean all
make CCC=${MAC_FULL_EXT} clean


# Create libmpi links for MAC
LIBMPI_FOR_MAC___=$(echo ${MACWORLD_MPI_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
echo ' '
echo 'MPI libs link for MAC =' ${LIBMPI_FOR_MAC___}
export LIBMPI_FOR_MAC___

# Create libblas links for MAC
LIBLAS_FOR_MAC___=$(echo ${MACWORLD_BLAS_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
echo 'Blas libs link for MAC =' ${LIBLAS_FOR_MAC___}
export LIBLAS_FOR_MAC___

# Create libatlas links for MAC
LIBATLAS_FOR_MAC___=$(echo ${MACWORLD_ATLAS_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
echo 'Atlas libs link for MAC =' ${LIBATLAS_FOR_MAC___}
export LIBATLAS_FOR_MAC___

# Create liblapack links for MAC
LIBLAPACK_FOR_MAC___=$(echo ${MACWORLD_LAPACK_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
echo 'Lapack libs link for MAC =' ${LIBLAPACK_FOR_MAC___}
export LIBLAPACK_FOR_MAC___

# Create libintel links for MAC
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "Intel" ]]
then
  LIBINTEL_FOR_MAC___=$(echo ${MACWORLD_INTEL_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
  echo 'Intel links for MAC =' ${LIBINTEL_FOR_MAC___}
else
  LIBINTEL_FOR_MAC___=""
  echo 'No Intel links for MAC'
fi  
export LIBINTEL_FOR_MAC___  

# Create libgfortran links for MAC
if [[ "${MACWORLD_SERCOMPIL_ENV}" == "GNU" ]]
then
  LIBGNU_FOR_MAC___=$(echo ${MACWORLD_GFORTRAN_LIBS} | sed 's%[^ ][^ ]*%-l&%g')
  echo 'GNU links for MAC =' ${LIBGNU_FOR_MAC___}
else
   LIBGNU_FOR_MAC___=""
  echo 'No GNU links for MAC' 
fi
export LIBGNU_FOR_MAC___

echo ' '


# Copy include files for use by MAC-based applications
make CCC=${MAC_FULL_EXT} cpinclude


# Compile in mac0 and mac2 modes
make CCC=${MAC_FULL_EXT} lib0
make CCC=${MAC_FULL_EXT} lib2


# Copy back the saved architecture-based installations
if [ -d ${__BACKUP_DIR} ]
  then
    for dir in ${__BACKUP_DIR}/*
      do
        if [ -d $dir ]
	  then       
	    cp -r ${dir} lib/
	  fi
      done
    rm -rf ${__BACKUP_DIR}
fi


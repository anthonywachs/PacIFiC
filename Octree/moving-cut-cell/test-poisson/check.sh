#/bin/bash

# Check test cases

if [ -z $1 ]; then
    for FILE in *.c; do
	export FILENAME="${FILE%.*}";

	#Check for no convergence
	grep -q "#" ${FILENAME}/log
	if [ $? == 0 ]; then
	    echo "Convergence problem in "${FILENAME}".c";
	fi
	
	# Check for difference with .ref
	if [ -f "${FILENAME}.ref" ];
	then
	    diff ${FILENAME}/log ${FILENAME}.ref
	    if [ $? -ne 0 ]; then
		echo ${FILENAME};
	    fi
	else
	    cp ${FILENAME}/log ${FILENAME}.ref
	fi
    done
else
    for FILE in $1.c; do
	export FILENAME="${FILE%.*}";
	
	#Check for no convergence
	grep -q "#" ${FILENAME}/log
	if [ $? == 0 ]; then
	    echo "Convergence problem in "${FILENAME}".c";
	fi
	
	# Check for difference with .ref
	if [ -f "${FILENAME}.ref" ];
	then
	    diff ${FILENAME}/log ${FILENAME}.ref
	    if [ $? -ne 0 ]; then
		echo ${FILENAME};
	    fi
	else
	    cp ${FILENAME}/log ${FILENAME}.ref
	fi
    done
fi

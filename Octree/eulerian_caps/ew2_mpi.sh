#!/bin/bash

# Step 1: update files on wachs2.math.ubc.ca
for header in /home/damien/phd/pacific/Octree/eulerian_caps/*.h;
do
	scp $header huet@wachs2.math.ubc.ca:/home/huet/Documents/phd/pacific/Octree/eulerian_caps/
done

for header in /home/damien/phd/pacific/Octree/eulerian_caps/navier-stokes/*.h;
do
        scp $header huet@wachs2.math.ubc.ca:/home/huet/Documents/phd/pacific/Octree/eulerian_caps/navier-stokes/
done

# for source_file in /home/damien/phd/pacific/Octree/eulerian_caps/devel-tests/*.c;
# do
#         scp $source_file huet@wachs2.math.ubc.ca:/home/huet/Documents/phd/pacific/Octree/eulerian_caps/devel-tests/
# done

for source_file in /home/damien/phd/pacific/Octree/eulerian_caps/devel-tests/*.c;
do
	s=$(basename $source_file)
	if [[ $2 == *"${s%.c}"* ]]; then
		scp $source_file huet@wachs2.math.ubc.ca:/home/huet/Documents/phd/pacific/Octree/eulerian_caps/devel-tests/$(basename $source_file)
	fi
done

scp /home/damien/phd/pacific/Octree/eulerian_caps/devel-tests/Makefile huet@wachs2.math.ubc.ca:/home/huet/Documents/phd/pacific/Octree/eulerian_caps/devel-tests/Makefile


# Step 2: launch simulation on wachs2.math.ubc.ca
ssh huet@wachs2.math.ubc.ca "cd Documents/phd/pacific ; source Env/PacIFiC-OpenMPI-1.10.2-GNU-4.8.5.env.sh > /dev/null 2>&1 ; cd Octree/eulerian_caps/devel-tests/ ; CC='mpicc -D_MPI=$1' make $2.tst"


# Step 3: copy results from wachs2.math.ubc.ca
scp -r huet@wachs2.math.ubc.ca:/home/huet/Documents/phd/pacific/Octree/eulerian_caps/devel-tests/$1 /home/damien/phd/pacific/Octree/eulerian_caps/devel-tests/$2

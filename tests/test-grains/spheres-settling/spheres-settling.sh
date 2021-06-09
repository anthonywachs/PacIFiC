#!/bin/bash

# 1) Run the test
./grains_seq Grains/Init/insert.xml 1> /dev/null 2> log

# 2) Define success
echo "success assessment:"

# 2.i) We first check that the maximum velocity is under 1.e-10 m/s
MAX_VEL=$(awk '{print $2}' Grains/Init/insert_VitesseMaxMean.dat | tail -1)
echo $(python3 -c $"if ${MAX_VEL}<1.e-10:"$'\n'$"    print(\"0\")"$'\n'$"else:"$'\n'$"    print(\"1\")") >> spheres-settling-success.txt
output="maximum velocity: "
output+=$(tail -n 1 spheres-settling-success.txt)
echo $output

# 2.ii) Then, we perform a systematic comparison of all output files
cd Reference_results
nb_assessments=1
for file in *
do
  let "nb_assessments++"
  cmp ../Grains/Init/$file $file ; echo $? >> ../spheres-settling-success.txt
  output="${file}: "
  output+=$(tail -n 1 ../spheres-settling-success.txt)
  echo $output
done
cd ..

# 2.iii) We run a command that fails if any of the above success assessments failed
[ $(grep -c "0" spheres-settling-success.txt) -eq $nb_assessments ]


set term qt size 1700,900

set xlabel "Angle"
set ylabel "x-component of the acceleration term"
set key samplen 2 spacing .8 font ",8" top left reverse invert Left
set xtics pi/4
set xrange [0:2*pi]
set format x '%.2Pπ'
set grid

isLayer1 = "((($1 < $8/2 + sqrt(2)*$7/2)) && ($1 > $8/2 - sqrt(2)*$7/2))"
isLayer2 = "((($1 < $8/2 + 3*sqrt(2)*$7/2)) && ($1 > $8/2 - 3*sqrt(2)*$7/2) && !(@isLayer1))"
isLayer3 = "((($1 < $8/2 + 5*sqrt(2)*$7/2)) && ($1 > $8/2 - 5*sqrt(2)*$7/2) && !(@isLayer2))"
isLayer4 = "((($1 > $8/2 + 5*sqrt(2)*$7/2)) || ($1 < $8/2 - 5*sqrt(2)*$7/2))"

set size 1,1
set origin 0,0
set multiplot layout 2,2 rowsfirst scale 0.9,0.9
set title "Layer 1: on the interface"
plot 'acceleration_1.csv' using (@isLayer1 ? (($2>0) ? $2 : $2+2*pi) : 1/0):4 with point lc -1 pt 6 title 'exact', \
'acceleration_1.csv' using (@isLayer1 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 1 title 'simulation - exact normals', \
'acceleration_2.csv' using (@isLayer1 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 2 title 'simulation - approx normals'

set title "Layer 2: 1 grid cell away from the interface"
plot 'acceleration_1.csv' using (@isLayer2 ? (($2>0) ? $2 : $2+2*pi) : 1/0):4 with point lc -1 pt 6 title 'exact', \
'acceleration_1.csv' using (@isLayer2 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 1 title 'simulation - exact normals', \
'acceleration_2.csv' using (@isLayer2 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 2 title 'simulation - approx normals'

set title "Layer 3: 2 grid cells away from the interface"
plot 'acceleration_1.csv' using (@isLayer3 ? (($2>0) ? $2 : $2+2*pi) : 1/0):4 with point lc -1 pt 6 title 'exact', \
'acceleration_1.csv' using (@isLayer3 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 1 title 'simulation - exact normals', \
'acceleration_2.csv' using (@isLayer3 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 2 title 'simulation - approx normals'

set title "Layer 4: 3 or more grid cells away from the interface"
plot 'acceleration_1.csv' using (@isLayer4 ? (($2>0) ? $2 : $2+2*pi) : 1/0):4 with point lc -1 pt 6 title 'exact', \
'acceleration_1.csv' using (@isLayer4 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 1 title 'simulation - exact normals', \
'acceleration_2.csv' using (@isLayer4 ? (($2>0) ? $2 : $2+2*pi) : 1/0):3 with point lc -1 pt 2 title 'simulation - approx normals'
unset multiplot


pause -1 "Hit any button to continue"

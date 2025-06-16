set xlabel 'Points per diameter'
set ylabel 'Error of the interface normals'
set title 'Convergence of interface normals in membrane cells for a circular capsule'
set logscale x
set logscale y
set format y '%h'
set format x '%h'
set grid
set key bottom left

set style line 1 lc -1 lw 1.5 dt 1 pt 4
set style line 2 lc -1 lw 1.5 dt 2 pt 5
set style line 3 lc -1 lw 1.5 dt 3 pt 6
set style line 4 lc -1 lw 1.5 dt 4 pt 7
set style line 5 lc -1 lw 1.5 dt 3 pt 8
set style line 6 lc -1 lw 1.5 dt 4 pt 9

plot 'out' using ($1 == 2 ? $2 : 1/0):($1 == 2 ? $3 : 1/0) w p ls 1 title "L_2-norm, height functions", \
'out' using ($1 == 2 ? $2 : 1/0):($1 == 2 ? $4 : 1/0) w p ls 2 title "L_{inf}-norm, height functions", \
'out' using ($1 == 1 ? $2 : 1/0):($1 == 1 ? $3 : 1/0) w p ls 3 title "L_2-norm, smoothed VOF", \
'out' using ($1 == 1 ? $2 : 1/0):($1 == 1 ? $4 : 1/0) w p ls 4 title "L_{inf}-norm, smoothed VOF", \
'out' using ($1 == 3 ? $2 : 1/0):($1 == 3 ? $3 : 1/0) w p ls 5 title "L_2-norm, exact level set", \
'out' using ($1 == 3 ? $2 : 1/0):($1 == 3 ? $4 : 1/0) w p ls 6 title "L_{inf}-norm, exact level set"


pause -1

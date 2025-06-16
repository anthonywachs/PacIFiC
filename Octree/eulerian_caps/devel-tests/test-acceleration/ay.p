set xlabel "Angle"
set ylabel "y-component of the acceleration term"
set key bottom right
set xtics pi/4
set format x '%.2Pπ'
set grid

plot 'acceleration.csv' using (($2>0) ? $2 : $2+2*pi):6 with point lc -1 pt 6 title 'exact', \
'acceleration.csv' using (($2>0) ? $2 : $2+2*pi):5 with point lc -1 pt 1 title 'simulation'

pause -1 "Hit any button to continue"

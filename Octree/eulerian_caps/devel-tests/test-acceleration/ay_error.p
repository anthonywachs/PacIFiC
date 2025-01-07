set xlabel "Angle"
set ylabel "Error on the y-component of the acceleration term (%)"
set xtics pi/4
set format x '%.2Pπ'
set grid

plot 'acceleration.csv' using (($2>0) ? $2 : $2+2*pi):((abs($6)>1.e-10) ? 100*($5-$6)/$6 : 1/0) with point lc -1 lw 1.2 pt 1 notitle

pause -1 "Hit any button to continue"

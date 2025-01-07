set y2tics
set ytics nomirror
set xlabel 'time'
set ylabel '|J-1|_{avg}'
set y2label '|J|_{max} - 1'
set key right
set grid
set title "J(t)"

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

plot 'out' every ::1 using 1:24 axis x1y1 w l ls 1 title 'lvl 6: |J-1|_{avg}', \
'out' every ::1 using 1:($8-1) axis x1y2 w l ls 2 title 'lvl 6: |J|_{max} - 1', \
'../sinusoidal-shear-lvl7/out' every ::1 using 1:24 axis x1y1 w l ls 3 title 'lvl7: |J-1|_{avg}', \
'../sinusoidal-shear-lvl7/out' every ::1 using 1:($8-1) axis x1y2 w l ls 4 title 'lvl 7: |J|_{max} - 1'

pause -1

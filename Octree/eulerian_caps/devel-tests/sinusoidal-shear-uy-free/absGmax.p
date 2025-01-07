set y2tics
set ytics nomirror
set xlabel 'time'
set ylabel '|G_{xx}|_{max} and |G_{yy}|_{max}'
set y2label '|G_{xy}|_{max}'
set key right
set grid
set title "Advection of G when u_x=y sin(2\pi t)"

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

plot 'out' every ::1 using 1:10 axis x1y1 w l ls 1 title '|G_{xx}|_{max}', \
'out' every ::1 u 1:14 axis x1y1 w l ls 2 title '|G_{yy}|_{max}', \
'out' every ::1 u 1:12 axis x1y2 w l ls 3 title '|G_{xy}|_{max}'

pause -1

set y2tics
set ytics nomirror
set xlabel 'time'
set ylabel '|T_s|_{max}'
set key right
set grid
set title "Stress tensor T_s when u_x=y sin(2\pi t)"

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

plot 'out' every ::1 using 1:16 axis x1y1 w l ls 1 title '|T_{s, xx}|_{max}', \
'out' every ::1 u 1:20 axis x1y1 w l ls 2 title '|T_{s, yy}|_{max}', \
'out' every ::1 u 1:18 axis x1y1 w l ls 3 title '|T_{s, xy}|_{max}'

pause -1

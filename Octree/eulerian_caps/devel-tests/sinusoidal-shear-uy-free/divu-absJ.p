set y2tics
set ytics nomirror
set xlabel 'time'
set ylabel '|J - 1|_{avg}'
set y2label '|div(u)|'
set key right
set grid
set title "Advection of J when u.y is not imposed"

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

back2 = back1 = 0
shift(x) = (back2 = back1, back1 = x)

plot 'out' every ::1 using 1:24 axis x1y1 w l ls 1 title '|J - 1|_{avg}', \
'out' every ::1 u 1:33 axis x1y2 w l ls 2 title 'div(u)_{avg}', \
'out' every ::1 u 1:32 axis x1y2 w l ls 3 title 'div(u)_{max}', \
'out' every ::1 u 1:36 axis x1y2 w l ls 4 title 'div(u_f)_{avg}', \
'out' every ::1 u 1:35 axis x1y2 w l ls 5 title 'div(u_f)_{max}'

pause -1

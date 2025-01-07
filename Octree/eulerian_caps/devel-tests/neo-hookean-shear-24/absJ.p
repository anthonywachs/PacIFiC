set xlabel 'time'
set ylabel '|J-1|_{avg}'
set grid

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

plot 'out' every ::1 using 1:24 w l ls 1 title 'G=t*t', \
'../neo-hookean-shear-25/out' every ::1 using 1:24 w l ls 2 title 'G advected', \
'../neo-hookean-shear-28/out' every ::1 using 1:24 w l ls 3 title 'G=t*t, ngcapsf in acceleration step', \
'../neo-hookean-shear-29/out' every ::1 using 1:24 w l ls 4 title 'G=t*t, ngcapsf in acceleration step, no extra mb layers', \

pause -1

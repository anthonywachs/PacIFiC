set xlabel 'time'
set ylabel '|J-1|_{avg}'
set ytics nomirror
set y2tics
set y2label '|a|_{avg}'
set grid

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

plot 'out' every ::1 using 1:24 axes x1y1 w l ls 1 title "|J-1|_{avg}", \
'out' every ::1 using 1:23 axes x1y2 w l ls 2 title "|a|_{avg}", \
'../neo-hookean-shear-41/out' every ::1 using 1:24 axes x1y1 w l ls 3 title "|J-1|_{avg} : NHS-41", \
'../neo-hookean-shear-41/out' every ::1 using 1:23 axes x1y2 w l ls 4 title "|a|_{avg} : NHS-41"

pause -1

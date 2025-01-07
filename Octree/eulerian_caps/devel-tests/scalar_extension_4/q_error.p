set xlabel "number of iterations"
set ylabel "L_2-error"
set y2label "L_{max}-error"
set y2tics nomirror

set logscale y
set logscale y2

set style line 1 lw 1.5 lc 'black' dt 1
set style line 2 lw 1.4 lc 'black' dt 2
set style line 3 lw 1.4 lc 'black' dt 3
set style line 4 lw 1.4 lc 'black' dt 4


plot 'out' every ::3 using 2:4 axis x1y1 w l ls 1 title "L_2-error, exact normals", \
'out' every ::3 using 2:5 axis x1y2 w l ls 2 title "L_{max}-error, exact normals", \
'../scalar_extension_5/out' every ::3 using 2:4 axis x1y1 w l ls 3 title "L_2-error, coarse normals", \
'../scalar_extension_5/out' every ::3 using 2:5 axis x1y2 w l ls 4 title "L_{max}-error, coarse normals"

pause -1

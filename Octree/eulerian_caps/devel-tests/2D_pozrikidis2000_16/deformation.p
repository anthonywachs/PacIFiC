set xlabel 'time'
set ylabel 'D'
set grid
set key bottom right

set style line 1 lc -1 lw 1.5 pt 4 pi 20
set style line 2 lc -1 lw 1.5 pt 6 pi 25 dt 2
set style line 3 lc -1 lw 1.5 pt 1 pi 30 dt 9
set style line 4 lc -1 lw 1.5 pt 12 pi 35 dt 8
set style line 5 lc -1 lw 1.5 pt 2 pi 28 dt 5

plot '/home/damien/phd/pacific/Octree/eulerian_caps/devel-tests/2D_pozrikidis2000/Ca0.003125.csv' using 1:($2/10) w l ls 1 title "Ca=0.003125 - Pozrikidis (2000)", \
'../2D_pozrikidis2000_16/out' using 1:23 w l ls 3 title "Ca=0.003125 - This study"


pause -1

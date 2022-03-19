
f(x) = 7./3.

plot f(x) w l lc rgb 'black' title 'theoretical', \
'stretch.txt' using 1 w l lc rgb 'black' dt 2 title 'computed - average stretch', \
'stretch.txt' using 2 w l lc rgb 'black' dt 3 title 'computed - minimum stretch', \
'stretch.txt' using 3 w l lc rgb 'black' dt 4 title 'computed - maximum stretch', \


pause -1

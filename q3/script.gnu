set palette gray
unset border
unset xtics
unset ytics
set pm3d
set xrange[0:100]
set yrange[0:100]
plot 'ex3.dat' using 1:2:3 w p lw 7 palette t ''

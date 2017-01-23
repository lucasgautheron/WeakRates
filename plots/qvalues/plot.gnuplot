set terminal epslatex color
set output 'plot.tex'

set pointsize 0.5

set xlabel '$Q_{exp}$ [MeV]'
set ylabel '$Q-Q_{exp}$ [MeV]'

set xrange [-20:20]
set yrange [-7.5:7.5]

plot '../../output/compare_qvalues.res' u ($9!=0?$9:NaN):($12-$9) w p ls 1 notitle, '' u ($9!=0?$9:NaN):($11-$9) w p ls 2 notitle, '' u ($9!=0?$9:NaN):($10-$9) w p ls 3 notitle, 1/0 w p ls 1 ps 1.5 t 'SEMF', 1/0 w p ls 2 ps 1.5 t 'DZ10', 1/0 w p ls 3 ps 1.5 t 'DZ33'
set output

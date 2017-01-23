set terminal epslatex color
set output 'plot.tex'

set pointsize 0.5

set xlabel 'Q [MeV]'
set ylabel '$\Delta Q$ [MeV]'

set yrange [-10:10]

plot '../../output/compare_qvalues.res' u ($9!=0?$9:NaN):($12-$9) t 'SEMF', '' u ($9!=0?$9:NaN):($11-$9) t 'DZ10', '' u ($9!=0?$9:NaN):($10-$9) t 'DZ33'
set output

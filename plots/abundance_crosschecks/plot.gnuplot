set terminal epslatex
set output 'plot.tex'

# set format xy "10^{%T}"
# set logscale xy

set xlabel '$\langle A \rangle$, $\langle Z \rangle$ EOS'
set ylabel '$\langle A \rangle$, $\langle Z \rangle$ Calc.'

set xrange [0:150]
set yrange [0:150]

set pointsize 0.2

plot '../../output/compare_neutrinos.res' u 5:7 every 100 lc rgb 'green' t 'A', '' u 6:8 every 100 lc rgb 'blue' t 'Z', x w l lc rgb 'black' notitle
set output

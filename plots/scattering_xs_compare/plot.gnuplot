set terminal epslatex
set output 'plot.tex'

set format xy "10^{%T}"
set logscale xy

set xrange [1e-48:1e-38]
set yrange [1e-48:1e-38]

set xlabel '$\sigma_{scatt}$ SNA [ cm$^{-2}$ ]'
set ylabel '$\sigma_{scatt}$ nuclei-distribution [ cm$^{-2}$ ]'

set pointsize 0.2

plot '../../output/compare_neutrinos.res' u 9:10 every 100 lc rgb 'red' t '', x w l lc rgb 'black' notitle
set output

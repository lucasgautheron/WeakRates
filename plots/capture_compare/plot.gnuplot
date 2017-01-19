set terminal epslatex
set output 'plot.tex'

set format xy "10^{%T}"
set logscale xy

set xrange [1e-40:1e-1]
set yrange [1e-40:1e-1]

set xlabel '$\lambda_{\mbox{ec}}$ SNA [ s$^{-1}$ ]'
set ylabel '$\lambda_{\mbox{ec}}$ distribution [ s$^{-1}$ ]'

set pointsize 0.05

plot '../../output/compare_capture.res' u 5:($1<100?$6:1/0) every 100 lc rgb 'red' t '', x w l lc rgb 'black' notitle
set output

set terminal epslatex color
set output 'plot.tex'

set format xy "10^{%T}"
set logscale xy

set xrange [1e-40:1e-1]
set yrange [1e-40:1e-1]

set xlabel '$\lambda_{\mbox{ec}}$ SNA [ s$^{-1}$ ]'
set ylabel '$\lambda_{\mbox{ec}}$ distribution [ s$^{-1}$ ]'

set pointsize 0.05

set palette model HSV rgbformulae 3,2,2

plot '../../output/compare_capture.res' u 5:($1<100?$6:1/0):2 every 100 t '' palette, x w l lc rgb 'black' notitle, x*100 w l dashtype 4 lc rgb 'black' notitle, x/100 w l dashtype 4 lc rgb 'black' notitle
set output

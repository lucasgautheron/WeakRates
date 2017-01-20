set terminal epslatex color
set output 'plot.tex'

set xlabel '$A$'
set ylabel '$Z$'

set pm3d map
set palette rgb 33,13,10

set zrange [0:9]

splot '../../output/compare_qvalues.res' u 1:2:(($3-$4)/$1) t 'B/A [MeV]'
set output

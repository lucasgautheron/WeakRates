ROOT=`pwd`
for d in `find . ! -path . -type d`
do
    cd "$ROOT/$d"
    echo `pwd`
    cp ../main.tex main.tex
    gnuplot plot.gnuplot
    pdflatex -interaction=batchmode main.tex
done

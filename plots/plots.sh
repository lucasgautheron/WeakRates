for d in `find . -type d`
do
    cd $d
    echo `pwd`
    gnuplot plot.gnuplot
    pdflatex -interaction=batchmode main.tex
done

# Set linestyle 1 to blue (#0060ad)
reset
set terminal
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.2

set title "Assignment unknown" font ",20"
set xtics font ", 20"
set ytics font ", 20"
set key font ",20"

plot 'STATIS' using 1:2 with linespoints linestyle 1  title "Exponential function"

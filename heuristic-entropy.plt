set terminal pdf
set output 'heuristic-entropy.pdf'

set xlabel "A"
set ylabel "S/N"

# Draw the critical point
set label at 9.87, 0.470311, 0 "" point pointtype 7 pointsize 0.25
set label at 9.87, 0.75, 0 "A_c"

plot './heuristic-entropy.dat' with lines

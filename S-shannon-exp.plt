set terminal pdf
set output 'S-shannon-exp.pdf'

set xlabel "N"
set ylabel "S"

plot './exp/S_shan/Sshan_A_13.00.dat' with linespoints, './exp/S_sub/Ssub_A_13.00.dat' with linespoints

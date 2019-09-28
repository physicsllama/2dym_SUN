set terminal pdf
set output 'S-shannon-fixedN.pdf'

set xlabel "N"
set ylabel "S"

plot './exp/S_shan/Sshan_N_21.dat' with lines, './exp/S_sub/Ssub_N_21.dat' with lines

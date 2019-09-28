set terminal pdf
set output 'S-shannon-brute.pdf'

set xlabel "N"
set ylabel "S"

plot './exp/S_shan/Sshan_A_10.00.dat' with lines, './exp/S_brute/Sbrute_A_10.00.dat' with lines

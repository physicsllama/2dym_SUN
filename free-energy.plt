set terminal pdf
set output 'free-energy.pdf'

set xlabel "A"
set ylabel "F"

plot './exp/F/F_DK.dat' with lines, './exp/F/F_DTV.dat' with lines, './exp/F/F_exp.dat' with lines, './exp/F/F_GM.dat' with lines, './exp/F/F_GS.dat' with lines
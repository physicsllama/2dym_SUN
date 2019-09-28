set terminal pdf
set output 'S-boltzmann.pdf'

set xlabel "A"
set ylabel "S"

f(x) = a*x**2 + b*x + c
fit f(x) "./exp/S_boltz_quad/Sboltz_A_10.00.dat" via a,b,c

plot "./exp/S_boltz_lin/Sboltz_A_10.00.dat", "./exp/S_boltz_quad/Sboltz_A_10.00.dat", f(x), b*x + c

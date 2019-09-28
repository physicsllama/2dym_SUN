set terminal pdf nocrop enhanced font "garamond,15"
set output 'rho.pdf'
set xlabel "position ({/Helvetica-Italic h})"
show xlabel
set ylabel "probability density"
show ylabel
set key

plot './exp/fermions/A5.00_N101.dat' title "A = 5.00", './exp/models/ratio_2.50.dat' with lines title "", './exp/fermions/A9.50_N101.dat' title "A = 9.50", './exp/models/ratio_4.75.dat' with lines title "", './exp/fermions/A15.00_N101.dat' title "A = 15.0", './exp/models/ratio_7.50.dat' with lines title ""

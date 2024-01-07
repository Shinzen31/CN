# Set output parameters
set terminal png
set output '/app/img/plot_RESVECGS.png'

# Set plot title and labels
set title "Convergence Plot for RESVECGS.dat"
set xlabel "Iteration"
set ylabel "Residual"

# Plot the data
plot '/app/data/RESVECGS.dat' with lines title 'RESVECGS'

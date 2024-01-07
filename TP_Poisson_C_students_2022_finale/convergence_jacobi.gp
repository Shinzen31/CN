# Set output parameters
set terminal png
set output '/app/img/plot_RESVECJac.png'

# Set plot title and labels
set title "Convergence Plot for RESVECJac.dat"
set xlabel "Iteration"
set ylabel "Residual"

# Plot the data
plot '/app/data/RESVECJac.dat' with lines title 'RESVECJac'

# Set output parameters
set terminal png
set output '/app/img/plot_RESVECalpha.png'

# Set plot title and labels
set title "Convergence Plot for RESVECalpha.dat"
set xlabel "Iteration"
set ylabel "Residual"

# Plot the data
plot '/app/data/RESVECalpha.dat' with lines title 'RESVECalpha'

using Plots
plotlyjs()  # Set the backend to PlotlyJS

x1 = LinRange(-8,8,101); x2 = LinRange(-8,8,101);
# f1(x1,x2) = 2 * x1^2 + 2 * x2^2 - 2 * x1 * x2;
f2(x1,x2) = (x1^4 - x2^4) * exp(0.1*(-x1^2-x2^2))
p = plot(x1,x2,f2,st=:surface,camera=(-60,40))

# plot(x1,x2,f,st=:contour,camera=(-30,60))
display(p)
# Save the plot to a file
savefig(p, "plots/static_plot.svg")
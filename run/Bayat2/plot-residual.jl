using Plots
using DelimitedFiles 
using LaTeXStrings

#pyplot()
plotlyjs() 
#pgfplotsx()
#gr()

# normal plot 
f= readdlm("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/Bayat2/residual0.txt")
#plot(f[:,1],f[:,3], label = ["force_x"])
#plot(f[:,2],f[:,4],label = ["force_y"])
plot(f[:,1],yaxis=:log) 
ylims!((1e-6,10000))


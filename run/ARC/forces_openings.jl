using Plots
using DelimitedFiles 
using LaTeXStrings

plotlyjs() 


ff= readdlm("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/ARC/forces0.txt")

plot(sqrt.(ff[:,1].*ff[:,1] +ff[:,2].*ff[:,2]),sqrt.(ff[:,3].*ff[:,3]+ff[:,4].*ff[:,4]),title = "t> to solution, tol xxx, dsip+=0.00xx, \n mesh X, XXX doffs" ,name = "reaction") 
plot!(sqrt.(ff[:,1].*ff[:,1] +ff[:,2].*ff[:,2]),sqrt.(ff[:,3].*ff[:,3]+ff[:,4].*ff[:,4]), seriestype=:scatter,markersize=:1,showlegend= false) 


plot!(sqrt.(ff[:,5].*ff[:,5] +ff[:,6].*ff[:,6]),sqrt.(ff[:,7].*ff[:,7]+ff[:,8].*ff[:,8]),name = "Mid plane magnitude") 

plot!(abs.(ff[:,5]),ff[:,7],name = "Mid plane x jump vs x force ") 
plot!(abs.(ff[:,6]),ff[:,8],name = "Mid plane y jump vs y force") 



plot!(legend=(0.77, +0.3))

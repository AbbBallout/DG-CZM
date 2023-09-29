using Plots
using DelimitedFiles 
using LaTeXStrings

plotlyjs() 


ff= readdlm("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/residual0.txt")
plot(ff[:,1], yaxis=:log , name = "residual0.txt") 

plot!(legend=(0.77, +0.7 ))

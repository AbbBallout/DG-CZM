using Plots
using DelimitedFiles 
using LaTeXStrings

plotlyjs() 


filelist = readdir("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs/q_points/") 
plot()
for file in filelist
    str="/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs/q_points/"*file
    local ff= readdlm(str)
    plot!(ff[:,3] ,ff[:,4],title = "At Q points" ,name = repr(ff[1,1]), mode="lines") 
    #plot!(ff[:,3] ,ff[:,5],title = "At Q points" ,name = repr(ff[1,1]), mode="lines") 
    

  #  plot!(abs.(ff[:,2]),ff[:,4],name = repr(ff[1,1])*" x jump vs x force") 
  #  plot!(abs.(ff[:,3]),ff[:,5],name = repr(ff[1,1])*" y jump vs y force") 
end



plot!(legend=(0.85, 1))




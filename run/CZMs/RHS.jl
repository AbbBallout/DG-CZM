using Plots
using DelimitedFiles 
using LaTeXStrings

plotlyjs() 


filelist = readdir("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs/system/") 
plot()
for file in filelist
     str="/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs/system/"*file
     if(occursin("rhs",str) )
      local ff= readdlm(str)
      #println(ff[1,3:end])
      plot!(ff[1,3:end]) 
     end
end

plot!(legend=(0.1, 1))




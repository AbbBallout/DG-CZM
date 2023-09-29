using Plots
using DelimitedFiles 
using LaTeXStrings

plotlyjs() 


filelist = readdir("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/q_points/") 
plot()
for file in filelist
    str="/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/q_points/"*file
    local ff= readdlm(str,skipstart=1)
    # g_eff_norm in the golabal system  vs  tr_norm() 
    plot!(ff[:,6] ,ff[:,9],title = "At Q points" ,name = repr(ff[1,2]), mode="lines") 
   # plot!(ff[:,4] ,ff[:,9] ,name = repr(ff[1,1]), mode="lines") 
   # plot!(ff[:,5] ,ff[:,9] ,name = repr(ff[1,1]), mode="lines") 
end

## PLoting the average of shear and normal stress
# ff= readdlm("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/q_points/q_point_1.txt",skipstart=1)
# tr_x=ff[:,7]
# tr_y = ff[:,8]
# mag=ff[:,9]
# counter=0
# for file in filelist
#     global counter=counter+1
#     str="/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/CZMs2/q_points/"*file
#     local ff= readdlm(str,skipstart=1)
#     global tr_x+=ff[:,7]
#     global tr_y+=ff[:,8]
#     global mag+=ff[:,9]

# end

# plot(ff[:,1],tr_x/counter,name="traction_x")
# plot!(ff[:,1],tr_y/counter,name="traction_y")
# plot!(ff[:,1],mag/counter,name="magnitude")


plot!(legend=(0.85, 0.7))




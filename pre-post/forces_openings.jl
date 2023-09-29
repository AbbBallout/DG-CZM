using Plots
using DelimitedFiles 
using LaTeXStrings

#pyplot()
plotlyjs() 
#pgfplotsx()
#gr()

# normal plot 
f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces1.txt")
plot(f[:,1],-f[:,2], label = ["force_x"])
plot!(f[:,1],f[:,3],label = ["force_y"])
plot!(f[:,1],sqrt.(f[:,2].*f[:,2]+f[:,3].*f[:,3]),title = "t> to solution, tol xxx, dsip+=0.00xx, \n mesh 2, XXX dofs" ,label = ["force_mag"]) 


#plot(f[:,1],f[:,4], label = ["interface_opening"])
#plot!(f[:,1],f[:,5], label = ["bulk_opening"])

## print the position of the max value 
#maxVal=findmax(f[:,2])
#println("maximum = ", maxVal[1]," achived at ",f[maxVal[2],1])

## print slopes 
#println((-f[2,2]+f[1,2])/(f[2,1]-f[1,1]))
#println((-f[end,2]+f[end-1,2])/(f[end,1]-f[end-1,1]))



## Plotting meshes 1 2 3 
# f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces0.txt")
# plot(f[:,1],sqrt.(f[:,2].*f[:,2]+f[:,3].*f[:,3]),label = "force - globaal.refine(4), ") 
# f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces1.txt")
# plot!(f[:,1],sqrt.(f[:,2].*f[:,2]+f[:,3].*f[:,3]) ,label = "force - globaal.refine(5), ")
# f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces2.txt") 
# plot!(f[:,1],sqrt.(f[:,2].*f[:,2]+f[:,3].*f[:,3]),title = "tol 1e-5,  dsip+=1e-4., \n cycles XXX " ,label = "force - globaal.refine(6) ") 


## interface opening 
# f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces0.txt")
# plot(f[:,1],f[:,4], label = ["interface_opening1"])
# plot!(f[:,1],f[:,5], label = ["bulk_opening1"])
# f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces1.txt")
# plot!(f[:,1],f[:,4], label = ["interface_opening2"])
# plot!(f[:,1],f[:,5], label = ["bulk_opening2"])
# f= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/forces2.txt")
# plot!(f[:,1],f[:,4], label = ["interface_opening3"])
# plot!(f[:,1],f[:,5], label = ["bulk_opening3"])
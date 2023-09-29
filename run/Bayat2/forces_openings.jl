using Plots
using DelimitedFiles 
using LaTeXStrings

#pyplot()
plotlyjs() 
#pgfplotsx()
#gr()

# normal plot 
ff= readdlm("/media/abbas/Ubuntu1/MainUbu/PhD-Parma/Solvers/Phase-field-Lornenzo-extended/run/Bayat2/forces0.txt")
#plot(ff[:,1],ff[:,3], label = ["fforce_x"])
#plot(ff[:,2],ff[:,4],label = ["fforce_y"])
plot(sqrt.(ff[:,1].*ff[:,1] +ff[:,2].*ff[:,2]),sqrt.(ff[:,3].*ff[:,3]+ff[:,4].*ff[:,4]),title = "t> to solution, tol xxx, dsip+=0.00xx, \n mesh X, XXX doffs" ,label = ["force_mag"]) 
plot!(sqrt.(ff[:,5].*ff[:,5] +ff[:,6].*ff[:,6]),sqrt.(ff[:,7].*ff[:,7]+ff[:,8].*ff[:,8])) 

plot!(sqrt.(ff[:,1].*ff[:,1] +ff[:,2].*ff[:,2]),sqrt.(ff[:,3].*ff[:,3]+ff[:,4].*ff[:,4]), seriestype=:scatter,markersize=:1 ) 

#plot(ff[:,1],ff[:,4], label = ["interfface_opening"])
#plot!(ff[:,1],ff[:,5], la
#label = ["bulk_opening"]

## print the position off the max value 
#maxVal=ffindmax(ff[:,2])
#println("maximum = ", maxVal[1]," achived at ",ff[maxVal[2],1])

## print slopes 
#println((-ff[2,2]+ff[1,2])/(ff[2,1]-ff[1,1]))
#println((-ff[end,2]+ff[end-1,2])/(ff[end,1]-ff[end-1,1]))



## Plotting meshes 1 2 3 
# ff= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-ffield-Lornenzo-extended/run/fforces0.txt")
# plot(ff[:,1],sqrt.(ff[:,2].*ff[:,2]+ff[:,3].*ff[:,3]),label = "fforce - globaal.reffine(4), ") 
# ff= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-ffield-Lornenzo-extended/run/fforces1.txt")
# plot!(ff[:,1],sqrt.(ff[:,2].*ff[:,2]+ff[:,3].*ff[:,3]) ,label = "fforce - globaal.reffine(5), ")
# ff= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-ffield-Lornenzo-extended/run/fforces2.txt") 
# plot!(ff[:,1],sqrt.(ff[:,2].*ff[:,2]+ff[:,3].*ff[:,3]),title = "tol 1e-5,  dsip+=1e-4., \n cycles XXX " ,label = "fforce - globaal.reffine(6) ") 


## interfface opening 
# ff= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-ffield-Lornenzo-extended/run/fforces0.txt")
# plot(ff[:,1],ff[:,4], label = ["interfface_opening1"])
# plot!(ff[:,1],ff[:,5], label = ["bulk_opening1"])
# ff= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-ffield-Lornenzo-extended/run/fforces1.txt")
# plot!(ff[:,1],ff[:,4], label = ["interfface_opening2"])
# plot!(ff[:,1],ff[:,5], label = ["bulk_opening2"])
# ff= readdlm("/media/abbas/Ubuntu/MainUbu/PhD-Parma/Solvers/Phase-ffield-Lornenzo-extended/run/fforces2.txt")
# plot!(ff[:,1],ff[:,4], label = ["interfface_opening3"])
# plot!(ff[:,1],ff[:,5], label = ["bulk_opening3"])
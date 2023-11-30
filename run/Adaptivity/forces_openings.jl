using Plots
using DelimitedFiles
using LaTeXStrings

plotlyjs()


function area_under_curve(x, f)
    n = length(x)
    total_area = 0.0

    for i in 2:n
        dx = x[i] - x[i-1]
        total_area += 0.5 * (f[i-1] + f[i]) * dx
    end

    return total_area
end

plot()


function plot_data(file_name, plot_1, plot_2, plot_3, max, in,style)
    ff = readdlm("main/run/Adaptivity/outputEN_10/$file_name.txt")
    # println("Area under reaction= ", area_under_curve(sqrt.(ff[:, 1] .* ff[:, 1] + ff[:, 2] .* ff[:, 2]), sqrt.(ff[:, 3] .* ff[:, 3] + ff[:, 4] .* ff[:, 4])))
    # println("Area under shear= ", area_under_curve(abs.(ff[:, 6]), ff[:, 7]))
    # println("Area under normal= ", area_under_curve(abs.(ff[:, 5]), ff[:, 8]))
    # println("Area under magnitude= ", area_under_curve(sqrt.(ff[:, 5] .* ff[:, 5] + ff[:, 6] .* ff[:, 6]), sqrt.(ff[:, 7] .* ff[:, 7] + ff[:, 8] .* ff[:, 8])))

    if (plot_1)
        x = sqrt.(ff[1:end, 1] .* ff[1:end, 1] + ff[1:end, 2] .* ff[1:end , 2])
        y = sqrt.(ff[1:end, 3] .* ff[1:end, 3] + ff[1:end, 4] .* ff[1:end, 4])
        plot!(x,
            y,
            #  color="blue"  ,
            linestyle=style,            
            linewidth=3,
            name=in)

        plot!(x, y,
            seriestype=:scatter,
            color="black", markersize=:1.5,
            showlegend=false
        )
    end

    if (plot_2)
        x = sqrt.(ff[1:end, 5] .* ff[1:end, 5] + ff[1:end, 6] .* ff[1:end, 6])
        y = sqrt.(ff[1:end, 7] .* ff[1:end, 7] + ff[1:end, 8] .* ff[1:end, 8])
        plot!(x, y,
            linewidth=1.5,
            color="red",
            name="Mid-plane force vs average jump (magnitude) "
        )

        # plot!(x, y,
        #     seriestype=:scatter,
        #     color="black",
        #     markersize=:1.5,
        #     showlegend=false)

    end

    if (plot_3)
        plot!(abs.(ff[1:end, 6]), ff[1:end, 7],
            #color="orange"  ,
            linewidth=1.5,
            name="Mid-plane shear force vs tangential jump")



        plot!(abs.(ff[:, 5]), ff[:, 8],
            #color="cyan"  ,
            linewidth=1.5,
            name="Mid-plane normal force vs normal jump")
    end
end

#plot_data("forces0", true, false, false, 0, "compress 0 ",:solid)
#plot_data("forces1", true, false, false, 0, "compress 0.1 ",:solid)
#plot_data("forces2", true, false, false, 0, "compress 0.2 ",:solid)
#plot_data("forces3", true, false, false, 0, "compress 0.3 ",:solid)
#plot_data("forces4", true, false, false, 0, "compress 0.4 ",:solid)
#plot_data("forces5", true, false, false, 0, "compress 0.5 ",:solid)
#plot_data("forces6", true, false, false, 0, "compress 0.6 ",:solid)
#plot_data("forces7", true, false, false, 0, "compress 0.7 ",:solid)




### Matrix fiber  

# plot_data("forces4", true, false, false, 162, "n_elements=4,908,  degree 1 -> dofs=39,264",:dash)
# plot_data("forces5", true, false, false, 162, "n_elements=20,064, degree 1 -> dofs=160,512",:dash)
# plot_data("forces6", true, false, false, 162, "n_elements=80,650, degree 1 -> dofs=645,200",:dash)

# plot_data("forces1", true, false, false, 162, "n_elements=4,908,  degree 2 -> dofs=88,344",:dot)
# plot_data("forces2", true, false, false, 162, "n_elements=20,064, degree 2 -> dofs=361,152",:dot)


#plot_data("mu=0", true, false, false, 0, "n_elements=4,908, degree 2 -> dofs=88,344 / eta=0",:dot)
# plot_data("mesh4degree1eta0.05", true, false, false, 0, "n_elements=80,650, degree 1 -> dofs=645,200 / eta=0",:solid)
# plot_data("mu=0.02", true, false, false, 250, "n_elements=20,064, degree 2 -> dofs=88,344 / eta=0.02",:dot)



### EN_10
plot_data("forces0", true, false, false, 0, " Base ",:solid)



#plot!(legend=(0.8, 0.4))
plot!(legend=:bottomright)
xlabel!("Displacement [mm]")
ylabel!("Force [N/mm]")
plot!(grid=true, gridlinewidth=3)
#plot!(size = [600,100])
#savefig("compression_plot.pdf")
#title!("reg3 made no difference")


# using BSplineKit

# Nguyenx = [0.0, 1.4e-3, 4.9e-3, 6e-3, 8.15e-3, 8.4e-3]
# Nguyeny = [0.0, 8.8, 18.9, 21.0, 17.8, 16.6]

# S = interpolate(Nguyenx, Nguyeny, BSplineOrder(4))
# xlims!(0, 8.15e-3)
# plot!(Nguyenx -> S(Nguyenx),
#     #  color="blue"  ,
#     linewidth=2.5,
#     seriestype=:line,
#     linestyle=:solid,
#     name="Nguyen V.P. 2014")
# xlims!(0, 0.0095)

# plot!(legend=(0.4, 0.6))
# plot!(size = [600,100])
# savefig("MF_plot_2.pdf")

######## ENF analytical ####
# L_enf= 70.0
# a_enf= 35.0
# E_enf= 1e+5
# h_enf = 3.0
# b_enf=  1.0
# I_enf= b_enf*(h_enf^3)/12
# min=30
# max=70
# P_enf = range(0, max, length=100)
# x_ENF = P_enf*(2*(L_enf^3)+3*(a_enf^3))/(96*E_enf*I_enf)
# plot!(x_ENF,P_enf)

# G_IIc_enf = 1
# P_enf = range(min, max, length=100)
# x_ENF = (P_enf/(96*E_enf*I_enf)).*( 2*L_enf^3 .+ (64*G_IIc_enf*b_enf*E_enf*I_enf)^1.5 ./ (sqrt(3)*P_enf.^3) )
# plot!(x_ENF,P_enf)

# P_enf = range(min+1, max, length=100)
# x_ENF = (P_enf/(24*E_enf*I_enf)).*( 2*L_enf^3 .- (64*G_IIc_enf*b_enf*E_enf*I_enf)^1.5 ./ (4*sqrt(3)*P_enf.^3) )
# plot!(x_ENF,P_enf)

# P_enf = range(min, max, length=100)
# x_ENF = P_enf*L_enf^3/(12*E_enf*I_enf)
# plot!(x_ENF,P_enf)


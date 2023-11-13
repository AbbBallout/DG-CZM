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


function plot_data(file_name, plot_1, plot_2, plot_3, max,in)
    ff = readdlm("main/run/Adaptivity/output/$file_name.txt")
    # println("Area under reaction= ", area_under_curve(sqrt.(ff[:, 1] .* ff[:, 1] + ff[:, 2] .* ff[:, 2]), sqrt.(ff[:, 3] .* ff[:, 3] + ff[:, 4] .* ff[:, 4])))
    # println("Area under shear= ", area_under_curve(abs.(ff[:, 6]), ff[:, 7]))
    # println("Area under normal= ", area_under_curve(abs.(ff[:, 5]), ff[:, 8]))
    # println("Area under magnitude= ", area_under_curve(sqrt.(ff[:, 5] .* ff[:, 5] + ff[:, 6] .* ff[:, 6]), sqrt.(ff[:, 7] .* ff[:, 7] + ff[:, 8] .* ff[:, 8])))

    if (plot_1)
        x = sqrt.(ff[1:end-max, 1] .* ff[1:end-max, 1] + ff[1:end-max, 2] .* ff[1:end-max, 2])
        y = sqrt.(ff[1:end-max, 3] .* ff[1:end-max, 3] + ff[1:end-max, 4] .* ff[1:end-max, 4])
        plot!(x,
            y,
            #  color="blue"  , 
            linewidth=2.5,
            name=in)

        plot!(x, y,
            seriestype=:scatter,
            color="black", markersize=:2,
            showlegend=false
        )
    end

    if (plot_2)
        x = sqrt.(ff[1:end-max, 5] .* ff[1:end-max, 5] + ff[1:end-max, 6] .* ff[1:end-max, 6])
        y = sqrt.(ff[1:end-max, 7] .* ff[1:end-max, 7] + ff[1:end-max, 8] .* ff[1:end-max, 8])
        plot!(x, y,
            linewidth=1.5,
            color="red",
            name="Mid-plane force vs average jump (magnitude) "
        )

        plot!(x, y,
            seriestype=:scatter,
            color="black",
            markersize=:1.5,
            showlegend=false)

    end

    if (plot_3)
        plot!(abs.(ff[:, 6]), ff[:, 7],
            #color="orange"  , 
            linewidth=1.5,
            name="Mid-plane shear force vs tangential jump")



        plot!(abs.(ff[:, 5]), ff[:, 8],
            #color="cyan"  , 
            linewidth=1.5,
            name="Mid-plane normal force vs normal jump")
    end
end


#plot_data("forces3.1", true, false, false ,0," STH")
#plot_data("forces3.2", true, false, false ,0," STH but again")
plot_data("forces1", true, true, true ,0," STH")






# This is data I constructed based on the values present in the table on the overleaf
# plot_data("forces1", true, false, false ,0," dofs 9632 , k=1e+5, intrinsic, Note:")
# plot_data("forces2", true, false, false ,0," dofs 9632 , k=1e+5, extrinsic, Note:")
# plot_data("forces3", true, false, false ,0," dofs 9632 , k=1e+6, intrinsic, Note:")
# plot_data("forces4", true, false, false ,0," dofs 9632 , k=1e+6, extrinsic, Note: ")
# plot_data("forces5", true, false, false ,0," dofs 9632 , k=1e+7, intrinsic, Note:")
# plot_data("forces6", true, false, false ,0," dofs 9632 , k=1e+7, extrinsic, Note: very premature")

#plot_data("forces7", true, false, false ,0," dofs 39264 , k=1e+5, intrinsic,Note:")
#plot_data("forces8", true, false, false ,0," dofs 39264 , k=1e+5, extrinsic,Note:")
#plot_data("forces9", true, false, false ,0," dofs 39264 , k=1e+6, intrinsic,Note:")
#plot_data("forces10", true, false, false ,0," dofs 39264 , k=1e+6, extrinsic,Note:")
#plot_data("forces11", true, false, false ,0," dofs 39264 , k=1e+7, intrinsic,Note:")
#plot_data("forces12", true, false, false ,0," dofs 39264 , k=1e+7, extrinsic,Note:")

# plot_data("forces13",  true, false, false ,0," dofs 160512 , k=1e+5, intrinsic, Note:")
# plot_data("forces14", true, false, false ,0," dofs 160512 , k=1e+5, extrinsic, Note:")
# plot_data("forces15",  true, false, false ,0," dofs 160512 , k=1e+6, intrinsic, Note:")
# plot_data("forces16", true, false, false ,0," dofs 160512 , k=1e+6, extrinsic, Note:")
# plot_data("forces17", true, false, false ,0," dofs 160512 , k=1e+7, intrinsic, Note:")
# plot_data("forces18", true, false, false ,0," dofs 160512 , k=1e+7, extrinsic, Note:")









plot!(legend=:topright)
xlabel!("Displacement [mm]")
ylabel!("Force [N/mm]")
plot!(grid=true, gridlinewidth=3)
#title!("reg3 made no difference")


# using BSplineKit

# Nguyenx = [0.0, 1.4e-3, 4.9e-3, 6e-3, 8.15e-3,8.4e-3]
# Nguyeny = [0.0, 8.8, 18.9, 21.0, 17.8,16.6]

# S = interpolate(Nguyenx, Nguyeny, BSplineOrder(4))
# xlims!(0, 8.15e-3)
# plot!(Nguyenx -> S(Nguyenx),
#     #  color="blue"  , 
#     linewidth=2.5,
#     seriestype=:line,
#     linestyle=:dash,
#     name="Nguyen V.P. 2014")
# xlims!(0,0.012)


#plot!(size = [600,100])
#plot!(legend=(0.50, 0.98))   
#savefig("MF_plot.pdf")

######### ENF analytical ####
# L_enf=70.0
# a_enf=40.0
# E_enf=1.0e+5
# h_enf = 3
# b_enf= 5
# I_enf= b_enf*h_enf^3/12
# min=20
# max=50+1
# P_enf = range(0, max, length=100)
# x_ENF = P_enf*(2*L_enf^3+3*a_enf^3)/(96*E_enf*I_enf)
# plot!(x_ENF,P_enf)

# G_IIc_enf = 0.03
# P_enf = range(min, max, length=100)
# x_ENF = (P_enf/(96*E_enf*I_enf)).*( 2*L_enf^3 .+ (64*G_IIc_enf*b_enf*E_enf*I_enf)^1.5 ./ (sqrt(3)*P_enf.^3) ) 
# plot!(x_ENF,P_enf)

# P_enf = range(min, 30, length=100)
# x_ENF = P_enf*L_enf^3/(12*E_enf*I_enf)
# plot!(x_ENF,P_enf)


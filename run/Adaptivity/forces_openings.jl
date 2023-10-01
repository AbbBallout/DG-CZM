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

function plot_data(file_name, plot_1, plot_2, plot_3)

    ff = readdlm("main/run/Adaptivity/output1/$file_name.txt")
    # println("Area under reaction= ", area_under_curve(sqrt.(ff[:, 1] .* ff[:, 1] + ff[:, 2] .* ff[:, 2]), sqrt.(ff[:, 3] .* ff[:, 3] + ff[:, 4] .* ff[:, 4])))
    # println("Area under shear= ", area_under_curve(abs.(ff[:, 6]), ff[:, 7]))
    # println("Area under normal= ", area_under_curve(abs.(ff[:, 5]), ff[:, 8]))
    # println("Area under magnitude= ", area_under_curve(sqrt.(ff[:, 5] .* ff[:, 5] + ff[:, 6] .* ff[:, 6]), sqrt.(ff[:, 7] .* ff[:, 7] + ff[:, 8] .* ff[:, 8])))

    if (plot_1)
        plot!(sqrt.(ff[:, 1] .* ff[:, 1] + ff[:, 2] .* ff[:, 2]),
            sqrt.(ff[:, 3] .* ff[:, 3] + ff[:, 4] .* ff[:, 4]),
            #  color="blue"  , 
            linewidth=1,
            name="Reaction on the top plane vs displacement")

        plot!(sqrt.(ff[:, 1] .* ff[:, 1] + ff[:, 2] .* ff[:, 2]), sqrt.(ff[:, 3] .* ff[:, 3] + ff[:, 4] .* ff[:, 4]),
            seriestype=:scatter,
            color="black",
            markersize=:1.25,
            showlegend=false
        )
    end

    if (plot_2)
    plot!(sqrt.(ff[:, 5] .* ff[:, 5] + ff[:, 6] .* ff[:, 6]),
        sqrt.(ff[:, 7] .* ff[:, 7] + ff[:, 8] .* ff[:, 8]),
        linewidth=1,
        color="red",
        name="Average jump vs force in magnitude"
    )

    plot!(sqrt.(ff[:, 5] .* ff[:, 5] + ff[:, 6] .* ff[:, 6]), sqrt.(ff[:, 7] .* ff[:, 7] + ff[:, 8] .* ff[:, 8]),
        seriestype=:scatter,
        color="black",
        markersize=:1.5,
        showlegend=false)

    end

    if(plot_3)
    plot!(abs.(ff[:, 6]), ff[:, 7],
        #color="orange"  , 
        linewidth=1,
        name="Average shear jump vs shear force")

    plot!(abs.(ff[:, 5]), ff[:, 8],
        #color="cyan"  , 
        linewidth=1,
        name="Average normal jump vs normal force ")
    end

end

plot_data("forces1", true, false, false)
plot_data("forces2", true, false, false)

xlabel!("Distance [mm]")
ylabel!("Force [N/m]")
#title!("Force vs Distance")

#plot!(legend=(0.55, 0.92))
plot!(legend=:bottomright)
 

# using BSplineKit

# Nguyenx = [0.0, 1.4e-3, 4.9e-3, 6e-3, 8.15e-3]
# Nguyeny = [0.0, 8.8, 18.9, 21.0, 17.8]

# S = interpolate(Nguyenx, Nguyeny, BSplineOrder(5))
# xlims!(0, 8.15e-3)
# plot!(Nguyenx -> S(Nguyenx),
#     #  color="blue"  , 
#     linewidth=1.5,
#     name="Nguyen")
# xlims!(0, 17e-3)

#savefig("Matrix_Fiber_seperation_plot.pdf")

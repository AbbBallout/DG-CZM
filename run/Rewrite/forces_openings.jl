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
    ff = readdlm("main/run/Rewrite/output/$file_name.txt")
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
            name="Mid-plane shear force vs AVG tangential jump")



        plot!(abs.(ff[1:end, 5]), ff[1:end, 8],
            #color="cyan"  ,
            linewidth=1.5,
            name="Mid-plane normal force vs AVG normal jump")
    end
end

plot_data("forces0", true, true, true, 0, " Base ",:solid)


#plot!(legend=(0.8, 0.4))
plot!(legend=:topright)
xlabel!("Displacement [mm]")
ylabel!("Force [N/mm]")
plot!(grid=true, gridlinewidth=3)
#plot!(size = [600,100])
#savefig("compression_plot.pdf")
#title!("reg3 made no difference")





using Plots
plotlyjs()


function heaviside(x)
    if x < 0
        return 0.0
    else
        return 1.0
    end
end

#### Warning t_n values are based on freddi's Ku. You can't use TCZn  

function analytical(mu_f_value, tn_value, us_range)

    K = [sigma_cs/delta_cs 0; 0 sigma_cn/delta_cn]

    u_n = 0
    if (tn_value <= 0)
        u_n = tn_value / (K[2, 2] - theta)
    else
        
        u_n = (tn_value / K[2, 2])
        
        
    end

    if (tn_value <= 0)
        u_n = -u_n
    end
    
    un = u_n * ones(size(us_range))
    u = [us_range, un]
    u_norm = sqrt.(u[1] .^ 2 + u[2] .^ 2)
    I = u[2] ./ u_norm
    II = u[1] ./ u_norm
    delta_c = 1 ./ sqrt.(((K[1, 1]) * II / sigma_cs) .^ 2 + ((K[2, 2]) * I / sigma_cn) .^ 2)
    delta_u = 1 ./ (K[1, 1] * delta_c .* II .^ 2 / (2 * GII) + K[2, 2] * delta_c .* I .^ 2 / (2 * GI))
    Damage = delta_u .* (u_norm - delta_c) ./ (u_norm .* (delta_u - delta_c))


    g = zeros(size(us_range))
    p = [zeros(size(us_range)), zeros(size(us_range))]
    q = [zeros(size(us_range)), zeros(size(us_range))]
    r = [zeros(size(us_range)), zeros(size(us_range))]

    for i in 2:length(us_range)
        tau = [K[1, 1] * (u[1][i-1] - (p[1][i-1] + q[1][i-1])) + r[1][i-1], K[2, 2] * (u[2][i-1] - (p[2][i-1] + q[2][i-1])) + r[2][i-1]]
        Coulumb = mu_f_value * tau[2] + abs.(tau[1])
        if (Coulumb > 0)
            g[i] = g[i-1] + Coulumb * sign(tau[1]) / K[1, 1]
        end
        p[1][i] = heaviside(u[2][i]) * u[1][i] + heaviside(-u[2][i]) * g[i]
        q[2][i] = heaviside(u[2][i]) * u[2][i]
        r[2][i] = heaviside(-u[2][i]) * theta * u[2][i]

        u_norm[i] = sqrt(u[1][i]^2 + u[2][i]^2)
        I[i] = u[2][i] / u_norm[i]
        II[i] = u[1][i] / u_norm[i]

        delta_c[i] = 1 / sqrt(((K[1, 1]) * II[i] / sigma_cs)^2 + ((K[2, 2]) * I[i] / sigma_cn)^2)
        delta_u[i] = 1 / (K[1, 1] * delta_c[i] * II[i]^2 / (2 * GII) + K[2, 2] * delta_c[i] * I[i]^2 / (2 * GI))

        Damage[i] = delta_u[i] * (u_norm[i] - delta_c[i]) / (u_norm[i] * (delta_u[i] - delta_c[i]))
        if (Damage[i] < 0)
            Damage[i] = 0
        end

        if (Damage[i] > 1)
            Damage[i] = 1
        end
    end


    TCZ = [K[1, 1] * (u[1] - Damage .* (p[1] + q[1])) + r[1], K[2, 2] * (u[2] - Damage .* (p[2] + q[2])) + r[2]]

    return un, TCZ
end

function plot_2D(mu_f_value, tn_value, us_range, un_range)

    plot!(us_range, TCZ[1], linewidth=2, name="mu= $mu_f_value tn_value=$tn_value")
    #plot!(un, TCZ[2], linewidth=2, name="mu= $mu_f_value tn_value=$tn_value")

end



plot()

plot!(legend=:topright)
xlabel!("u_s [mm]")
ylabel!("sigma_s [MPa]")
plot!(grid=true, gridlinewidth=3)

us = 0:0.0001:1.2

sigma_cs = 1.0
delta_cs = 0.1
sigma_cn = 1.0
delta_cn = 0.1
GII = 0.5
GI = 0.5
theta = 1e+5

# un, TCZ = analytical(0,-1,us)
# plot_2D(0,-1,us,un)
# un, TCZ = analytical(0.2,-1,us)
# plot_2D(0.2,-1,us,un)
# un, TCZ = analytical(0.4,-1,us)
# plot_2D(0.4,-1,us,un)
# un, TCZ = analytical(0.6,-1,us)
# plot_2D(0.6,-1,us,un)
# un, TCZ = analytical(0.8,-1,us)
# plot_2D(0.8,-1,us,un)

#ylims!(0, 1.4)
#savefig("analytical_mu_plot.svg")

# un, TCZ = analytical(0.5,1.5,us)
# plot_2D(0.5,1.5,us,un)
# un, TCZ = analytical(0.5,0.75,us)
# plot_2D(0.5,0.75,us,un)
# un, TCZ = analytical(0.5,0,us)
# plot_2D(0.5,0,us,un)
# un, TCZ = analytical(0.5,-0.75,us)
# plot_2D(0.5,-0.75,us,un)
# un, TCZ = analytical(0.5,-1.5,us)
# plot_2D(0.5,-1.5,us,un)

# ylims!(0, 1.4)
#savefig("analytical_tn_value_plot.svg")



mu_f = 0.6
us = 0:0.01:1.25
t_n_range = -1:0.01:10
TCZ_3D_s = zeros(length(us), length(t_n_range))
TCZ_3D_n = zeros(length(us), length(t_n_range))
u_n_vals = zeros(length(t_n_range))

for i in 1:length(t_n_range)

    local un, TCZ = analytical(mu_f, t_n_range[i], us)
    TCZ_3D_s[:,i]= TCZ[1]
    TCZ_3D_n[:,i]= TCZ[2]
    u_n_vals[i] = un[1]

end

using LaTeXStrings
gr()
surface(u_n_vals, us, TCZ_3D_s, c = :jet)
plot!(
    xlabel = L"[\![u_n ]\!] \, mm",
    ylabel = L"[\![u_s ]\!] \, mm",
    zlabel = L"\tau_s \, MPa",
    #title = "",
    c = TCZ_3D_s,  # Use data values for color
    linetype = :auto,  # Set line type, e.g., :auto, :surface, :wireframe
    alpha = 0.1,  # Set line transparency
    grid = true,  # Disable grid lines
    legend = :outertopright,  # Move legend to the outer top-right corner
    camera = (70, 35),  # Adjust camera angles
)

# savefig("analytical_TCZ3D.svg")


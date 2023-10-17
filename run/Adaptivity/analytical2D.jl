using Plots
plotlyjs()


function heaviside(x)
    if x < 0
        return 0.0
    else
        return 1.0
    end
end


function analytical(mu_f, tn)

    us = 0:0.0001:1.2

    sigma_cs = 1.0
    delta_cs = 0.1
    sigma_cn = 1.0
    delta_cn = 0.1

    K = [sigma_cs/delta_cs 0; 0 sigma_cn/delta_cn]

    GII = 0.5
    GI = 0.5


    u_n = 0
    theta = 1e+4
    if (tn <= 0)
        u_n = tn / (K[2, 2] - theta)
    else
        delta_uu = 2 * GI / sigma_cn

        u_n = (tn / K[2, 2])
        if (u_n > delta_uu)
            u_n = delta_uu
        end
    end

    if (tn <= 0)
        u_n = -u_n
    end

    un = u_n * ones(size(us))
    u = [us, un]
    u_norm = sqrt.(u[1] .^ 2 + u[2] .^ 2)
    I = u[2] ./ u_norm
    II = u[1] ./ u_norm
    delta_c = 1 ./ sqrt.(((K[1, 1]) * II / sigma_cs) .^ 2 + ((K[2, 2]) * I / sigma_cn) .^ 2)
    delta_u = 1 ./ (K[1, 1] * delta_c .* II .^ 2 / (2 * GII) + K[2, 2] * delta_c .* I .^ 2 / (2 * GI))
    Damage = delta_u .* (u_norm - delta_c) ./ (u_norm .* (delta_u - delta_c))


    g = zeros(size(us))
    p = [zeros(size(us)), zeros(size(us))]
    q = [zeros(size(us)), zeros(size(us))]
    r = [zeros(size(us)), zeros(size(us))]

    for i in 2:length(us)
        tau = [K[1, 1] * (u[1][i-1] - (p[1][i-1] + q[1][i-1])) + r[1][i-1], K[2, 2] * (u[2][i-1] - (p[2][i-1] + q[2][i-1])) + r[2][i-1]]
        Coulumb = mu_f * tau[2] + abs.(tau[1])
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


    plot!(us, TCZ[1], linewidth=2, name="mu= $mu_f tn=$tn")

end



plot()

plot!(legend=:topright)
xlabel!("u_s [mm]")
ylabel!("sigma_s [MPa]")
plot!(grid=true, gridlinewidth=3)


# analytical(0,-1)
# analytical(0.2,-1)
# analytical(0.4,-1)
# analytical(0.6,-1)
# analytical(0.8,-1)

# ylims!(0, 1.4)

# savefig("analytical_mu_plot.svg")


# analytical(0.5,1.5)
# analytical(0.5,0.75)
# analytical(0.5,0)
# analytical(0.5,-0.75)
# analytical(0.5,-1.5)

# ylims!(0, 1.4)
# savefig("analytical_tn_plot.svg")








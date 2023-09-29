using Plots 
using SymPy
#plotlyjs()
gr() 


@vars gs gn t0 lambdaf m  real=True

g=[gs,gn]
lambda=g.norm()
TCZ=t0*g*((lambdaf-lambda)/lambdaf)^m/lambda

TCZ=TCZ.subs([(t0,3),(lambdaf,0.4),(m,1)])

#gs = range(0, stop=0.4, length=50)
#gn = range(0, stop=0.4, length=50)


r_max = 0.4
θ_vals = range(0, stop=π/2, length=10)


gs_vals = []
gn_vals = []
for θ in θ_vals
    for r in range(0, stop=r_max, length=10)
        gs = r * cos(θ)
        gn = r * sin(θ)
        push!(gs_vals, gs)
        push!(gn_vals, gs)
    end
end



surface(gs_vals, gn_vals, TCZ[1])

using Plots
using SymPy

@vars gs gn t0 lambdaf m  real=True

g=[gs,gn]
lambda=g.norm()
TCZ=t0*g*((lambdaf-lambda)/lambdaf)^m/lambda

#TCZ[1]=t0*gs*((lambdaf-gs)/lambdaf)^m/gs
#TCZ[2]=t0*gn*((lambdaf-gn)/lambdaf)^m/gn



TCZ=TCZ.subs([(t0,3),(lambdaf,0.4),(m,1)])
Jaccobian = [diff(TCZ[1],gs) diff(TCZ[1],gn) ; diff(TCZ[2],gs) diff(TCZ[2],gn) ]

u = [0.1,0.1] # initial guess for u
t = [1.37132,1.37132] # where you want to converge

iteration=0

println("Start coupled")
while (iteration<5)
    global iteration+=1
    invJacobian=Jaccobian.subs([(gs,u[1]),(gn,u[2])])
    invJacobian=inv.(invJacobian)
    println("invJac=",invJacobian)
    #invJacobian[1,2]=0
    #invJacobian[2,1]=0 
    res=TCZ.subs([(gs,u[1]),(gn,u[2])])-t
    println("res=",res)
    deltau=-invJacobian*res
    u=u+deltau
    println("u=",u)
    println()
end

println("Start staggered")
while (iteration<5)
    global iteration+=1
    invJacobian=Jaccobian.subs([(gs,u[1]),(gn,u[2])])
    invJacobian=inv.(invJacobian)
    println("invJac=",invJacobian)
    #invJacobian[1,2]=0
    #invJacobian[2,1]=0 
    res=TCZ.subs([(gs,u[1]),(gn,u[2])])-t
    println("res=",res)
    deltau=-invJacobian*res
    u=u+deltau
    println("u=",u)
    println()
end




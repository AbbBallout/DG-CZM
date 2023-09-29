using Plots

n=100
t0=3
lambda_f=0.4
m=1


x0=0
x1=0.4-1e-3

u_jump=LinRange(x0, x1, n+1)
t_analytical=t0*(1 .- u_jump./lambda_f).^m
plot(u_jump,t_analytical)


res=0
iteration=0
u=0    ##initial guess for u
t_1=2  ## where you want to converge


while(iteration<3)

    global iteration+=1
    tangent=-m*t0/lambda_f*(1 - u/lambda_f)^(m-1) 
    traction=t0*(1 - u/lambda_f)^m - t_1 
    update=-traction/tangent
    global u=u+update
    
    println(traction)
    println(u)
    println()
   
end 





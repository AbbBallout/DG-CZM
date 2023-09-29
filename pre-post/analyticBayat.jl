using SymPy
using Plots
plotly

x0=1e-5
x1= 0.4

gn = symbols("gn", real=True)
gs = symbols("gs", real=True)
beta=symbols("beta")
geffs=beta*gs
geffn=(abs(gn)+gn)/2
lambda=sqrt(geffn^2+(geffs)^2)
lambdaf=symbols("lambdaf")
t_0=symbols("t_0")
m=symbols("m")
lambda_s=t_0*geffs*((lambdaf-lambda)/lambdaf)^m/lambda
lambda_n=t_0*geffn*((lambdaf-lambda)/lambdaf)^m/lambda

lambda_n=lambda_n.subs(t_0,3);
lambda_n=lambda_n.subs(lambdaf,0.4)
lambda_n=lambda_n.subs(beta,1)
lambda_n=lambda_n.subs(gs,0)
lambda_n=lambda_n.subs(m,1)
plot(lambda_n(gn),x0,x1)
print((lambda_n.subs(gn,x1)-lambda_n.subs(gn,x0))/(x1-x0)  )

CZM11=diff(lambda_n,gn)
plot!(CZM11,x0,x1)

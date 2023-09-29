using SymPy
using Plots
plotlyjs() 

x0= 0.01
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
plot(lambda_n,x0,x1)

CZM=diff(lambda_n,gn)
CZM=CZM.subs(t_0,3);
CZM=CZM.subs(lambdaf,0.4)
CZM=CZM.subs(beta,1)
CZM=CZM.subs(gs,0)
CZM=CZM.subs(m,1)
plot!(CZM,x0,x1,labels="CZM_analytical dervative")

#### The CZM in the paper ##########
KCZ_pp = 0.5*(gn/2-abs(gn)/2)^2
KCZ_pp= -KCZ_pp*(lambdaf+lambda*(m-1))/(lambda^2*(lambdaf-lambda))
KCZ_pp= KCZ_pp + 
KCZ_pp=KCZ_pp*(t_0*((lambdaf-lambda)/lambdaf)^m)/lambda
 KCZ_pp=KCZ_pp.subs(t_0,3);
 KCZ_pp=KCZ_pp.subs(lambdaf,0.4)
 KCZ_pp=KCZ_pp.subs(beta,1)
 KCZ_pp=KCZ_pp.subs(gs,0)
 KCZ_pp=KCZ_pp.subs(m,1)
plot!(KCZ_pp,x0,x1,labels="CZM_paper")

#### The CZM corrected ##########
KCZ_ppc = 0.5*(gn/2+abs(gn)/2)^2
KCZ_ppc= -KCZ_ppc*(lambdaf+lambda*(m-1))/(lambda^2*(lambdaf-lambda))
 KCZ_ppc= KCZ_ppc + beta
 KCZ_ppc=KCZ_ppc*(t_0*((lambdaf-lambda)/lambdaf)^m)/lambda
 KCZ_ppc=KCZ_ppc.subs(t_0,3);
 KCZ_ppc=KCZ_ppc.subs(lambdaf,0.4)
 KCZ_ppc=KCZ_ppc.subs(beta,1)
 KCZ_ppc=KCZ_ppc.subs(gs,0)
 KCZ_ppc=KCZ_ppc.subs(m,1)
plot!(KCZ_ppc,x0,x1,labels="CZM_corrected")

######
#lambda=lambda.subs(gs,0)
#plot(lambda,x0,x1)
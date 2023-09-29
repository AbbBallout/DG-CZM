using Plots 
using SymPy

@vars t_0  g  real=True

lambda_0=1e-2
lambda_f=0.4

TCZ1=t_0*g/lambda_0 
TCZ2=t_0*(lambda_f-g)/(lambda_f-lambda_0)

blending=g/lambda_0

TCZ=TCZ1*(1-blending) + TCZ2*(blending) 

TCZ=TCZ.subs([(t_0,3)])

plot(TCZ,0,lambda_0)

TCZ2=TCZ2.subs([(t_0,3)])
plot!(TCZ2,lambda_0,lambda_f)
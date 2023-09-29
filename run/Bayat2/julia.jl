using SymPy
using Plots
plotlyjs() 

x0= 0.05
x1= 0.4

g_n = symbols("g_n", real=True)
g_s = symbols("g_s", real=True)
beta=symbols("beta")
geff_s=beta*g_s
geff_n=(abs(g_n)+g_n)/2
lambda=sqrt(geff_n^2+geff_s^2)
lambda_f=symbols("lambda_f")
t_0=symbols("t_0")
m=symbols("m")
t_s=t_0*geff_s*((lambda_f-lambda)/lambda_f)^m/lambda
t_n=t_0*geff_n*((lambda_f-lambda)/lambda_f)^m/lambda



#t_n=t_n.subs(t_0,3);
#t_n=t_n.subs(lambda_f,0.4)
#t_n=t_n.subs(beta,1)
#t_n=t_n.subs(g_s,0)
#t_n=t_n.subs(m,1)
#plot(t_n,x0,x1)

CZM=diff(t_s,g_s)
#plot!(CZM,x0,x1,labels="CZM_analytical dervative")

# expression=(lambda_f+lambda(m-1))/(lambda^2*(lambda_f-lambda))
# expression=-expression*((abs(g_n)+g_n)/2)^2
# expression=expression+1
# expression=expression*t_0*(lambda_f/lambda_f-lambda/lambda_f)^m/lambda
# expression=expression.subs(t_0,3);
# expression=expression.subs(m,1);
# expression=expression.subs(lambda_f,0.4);
# expression=expression.subs(g_s,0);
# plot(expression,x0,x1,labels="expression")

###
expression=(lambda_f+lambda(m-1))/(lambda^2*(lambda_f-lambda))
expression=-expression*beta^3*(g_s)^2
expression=expression+1
expression=expression*t_0*(lambda_f/lambda_f-lambda/lambda_f)^m/lambda
expression=expression.subs(t_0,3);
expression=expression.subs(m,1);
expression=expression.subs(lambda_f,0.4);
expression=expression.subs(g_n,0);
expression=expression.subs(beta,1);
plot(expression,x0,x1,labels="expression")
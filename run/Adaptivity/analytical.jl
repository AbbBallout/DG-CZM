using SymPy

@vars us un real=true
u = Sym[us ; un]

@vars theta  g  real=true
p = [Heaviside(un)*us+Heaviside(-un)*g;0]
q=  [0 ; un*Heaviside(un)]

@vars  sigma_cs delta_cs sigma_cn delta_cn real=true
K = Sym[sigma_cs/delta_cs  0 ; 0  sigma_cn/delta_cn]
tu= K*u

@vars k_us k_un   real=true
#mixity=u[1]/u.norm();
k_u=(abs(u[1])/(abs(u[1])+abs(u[2])))*k_us + (abs(u[2])/(abs(u[1])+abs(u[2])))*k_un
sigma_c=(abs(u[1])/(abs(u[1])+abs(u[2])))*sigma_cs + (abs(u[2])/(abs(u[1])+abs(u[2])))*sigma_cn
Damage = k_u*(sigma_c - tu.norm())/( tu.norm()*(sigma_c-k_u))

TCZ= K*(u-Damage*(p+q))

Jac = Sym[1  0 ; 0  1] # dTCZ/du
Jac+= kron(-(p+q), Sym[diff(Damage,us)  diff(Damage,un)] ) 
Jac+= -Damage *  Sym[diff(p[1],us) diff(p[1],un) ; diff(p[2],us) diff(p[2],un)] 
Jac+= -Damage *  Sym[diff(q[1],us) diff(q[1],un) ; diff(q[2],us) diff(q[2],un)] 
Jac=Jac*K
#@vars ku sigma_c delta_c  real=true
#Jac=Jac.subs( [(kus, ku), (kun, ku) ,(sigma_cs,sigma_c),(sigma_cn,sigma_c),(delta_cs, delta_c),(delta_cn, delta_c) ] )

@vars mu 
tau_trial = K*(u-(p+q))
g = (1 / K[1][1]) * (mu * tau_trial[2] + abs(tau_trial[1])) * (abs(tau_trial[1]) / tau_trial[1]);


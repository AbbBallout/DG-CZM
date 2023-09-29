using SymPy

@vars us un real=true
u = Sym[us ; un]
lambda= u.norm()

# @vars kt kn delta_c  lambda_f 
# d= lambda_f*(lambda-delta_c)/(lambda*(lambda_f-delta_c))
# TCZ=Sym[kt*(1-d)*u[1] ; kn*(1-d)*u[2]]

#@vars sig_c delta_c n lambda_f  
#TCZ=sig_c*((lambda_f - lambda)/(lambda_f-delta_c))^n * u/lambda

## Stress damage based 
@vars kt kn ku sigma_c penalty n

K = Sym[kt 0 ; 0 kn]
sigma_eff= K*u
sigma_eq = sigma_eff.norm()
damage = ku*(sigma_c -sigma_eq)/(sigma_eq*(sigma_c-ku))
p = Sym[us ; 0]
q = Sym[0;Heaviside(un)*un] 
r = Sym[0,penalty*(-un)^n*Heaviside(-un)]

TCZ= K*(u-damage*(p+q)) 

Jac=Sym[diff(TCZ[1],us) diff(TCZ[1],un) ; diff(TCZ[2],us) diff(TCZ[2],un)] 

Jac_p = K # dTCZ/dg
Jac_p+= kron(-K*(p+q), Sym[diff(damage,us)  diff(damage,un)] )     
Jac_p+= -K*damage *  Sym[diff(p[1],us) diff(p[1],un) ; diff(p[2],us) diff(p[2],un)] 
Jac_p+= -K*damage *  Sym[diff(q[1],us) diff(q[1],un) ; diff(q[2],us) diff(q[2],un)] 


Jac[1,1] 
#Jac_p[1,1]

#Jac[1,2].subs( [(us, 1.3), (un, 0.55) ,(kt,0.32),(kn,15),(sigma_c,4.33),(ku,2.12) ] )
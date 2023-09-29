using SymPy

@vars us un real=true
u = Sym[us ; un]
lambda= u.norm()

# @vars kt kn delta_c  lambda_f 
# d= lambda_f*(lambda-delta_c)/(lambda*(lambda_f-delta_c))
# TCZ=Sym[kt*(1-d)*u[1] ; kn*(1-d)*u[2]]

#@vars sig_c delta_c n lambda_f  
#TCZ=sig_c*((lambda_f - lambda)/(lambda_f-delta_c))^n * u/lambda

################ 
#Stress damage based impliict
############################
# @vars kt kn ku sigma_c penalty n mu gn real =true

# K = Sym[kt 0 ; 0 kn]
# sigma_eff= K*u
# sigma_eq = sigma_eff.norm()
# damage = ku*(sigma_c -sigma_eq)/(sigma_eq*(sigma_c-ku))
# #p = Sym[us ; 0] 
# #or
# p = Sym[mu*(kn/kt)*(un)*abs(us-gn)/(us-gn) + us-gn; 0]
# q = Sym[0;Heaviside(un)*un] 
# r = Sym[0,penalty*(-un)^n*Heaviside(-un)]

# TCZ= K*(u-damage*(p+q))

# Jac=Sym[diff(TCZ[1],us) diff(TCZ[1],un) ; diff(TCZ[2],us) diff(TCZ[2],un)] 

# Jac_p = K # dTCZ/dg
# Jac_p+= kron(-K*(p+q), Sym[diff(damage,us)  diff(damage,un)] )     
# Jac_p+= -K*damage *  Sym[diff(p[1],us) diff(p[1],un) ; diff(p[2],us) diff(p[2],un)] 
# Jac_p+= -K*damage *  Sym[diff(q[1],us) diff(q[1],un) ; diff(q[2],us) diff(q[2],un)] 


# Jac[1,1] 
#Jac_p[1,1]

#Jac[1,2].subs( [(us, 1.3), (un, 0.55) ,(kt,0.32),(kn,15),(sigma_c,4.33),(ku,2.12) ] )

#TCZ= TCZ.subs( [(un, 1e-3) ,(kt,100),(kn,100),(sigma_c,0.1),(ku,0.5) ] )
#Jac_p= Jac_p.subs( [(un, -1e-3) ,(kt,100),(kn,100),(sigma_c,0.1),(ku,0.5) ] )

#plot(TCZ[1],0,5e-3)

################ 
#Stress damage based explicit
############################

# @vars kt kn lambda_f real =true
# K = Sym[kt 0 ; 0 kn]
# damage = lambda/lambda_f
# TCZ=K*(1-damage)*u

# Jac=Sym[diff(TCZ[1],us) diff(TCZ[1],un) ; diff(TCZ[2],us) diff(TCZ[2],un)] 
# Jac[2,2] 

# @vars x
# using Plots
# damage = 10*(0.1 -x)/(x*(0.1-10))
# plot(damage,[0.1:0.1:10])
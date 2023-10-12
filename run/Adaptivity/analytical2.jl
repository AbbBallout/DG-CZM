using SymPy
using Plots
gr()

@vars us un real=true
u = Sym[us ; Heaviside(un)*un]
 
sigma_cs= 1.2
delta_cs= 0.2
sigma_cn= 1.2 
delta_cn= 0.2

K = [sigma_cs/delta_cs  0 ; 0  sigma_cn/delta_cn]

GII = 0.5
GI = 0.5

II=u[1]/u.norm()
I =u[2]/u.norm() 

delta_c = 1/sqrt(((K[1,1])*II/sigma_cs)^2 + ((K[2,2])*I/sigma_cn)^2 )
delta_u = 1/(K[1,1]*delta_c*II^2/(2*GII) + K[2,2]*delta_c*I^2/(2*GI) )

Damage = delta_u*(u.norm()-delta_c)/(u.norm()*(delta_u-delta_c))

#Damage = max(0,min(1,Damage))
theta = 0.0
p=Sym[u[1];0]
q=Sym[0;Heaviside(u[2])*un] 
r=Sym[0,Heaviside(-u[2])*theta*un]  

TCZ1 = K*u  +r
TCZ2 = K*(u-Damage*(p+q)) +r 

number = 0.12
TCZ1=TCZ1.subs(un,number)
TCZ2=TCZ2.subs(un,number)
Damage=Damage.subs(un,number)

max = 0.82
vec=1e-8:0.001:max
position = 0  
fail = 0 
for i in 1:length(vec)
    if(Damage.subs(us,vec[i])>0)
        global position=i
        global fail=vec[i]
        break
    end
end



plot(TCZ1[1],0,fail)
plot!(TCZ2[1],fail,max)


plot!([0,number],[0,TCZ1[2]])

un_plt = zeros(size(vec)[1]-position)
TCZn_plt = zeros(size(vec)[1]-position)
for i in 1:length(vec)-position
    un_plt[i] = number  
    TCZn_plt[i] = TCZ2[2].subs(us,vec[i+position])
end

 

plot!(un_plt,TCZn_plt)








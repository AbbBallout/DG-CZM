using SymPy
using Plots

@vars gn

epsi=1e-6
lambda_f=0.1
sig_c=1

tn= sig_c*(lambda_f-gn-epsi)/lambda_f * gn/(gn+epsi)
plot(tn,0,lambda_f,labels="expression")
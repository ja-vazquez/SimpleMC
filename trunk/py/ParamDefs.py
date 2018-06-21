##
## This class has parameter defintions for all
## parameter used in this code.
##
## Change here for bounds, or import and rewrite.
##
##


from Parameter import *

## Parameters are value, variation, bounds
#0.3038, 0.02234, 0.6821
Om_par = Parameter("Om", 0.3038, 0.05, (0.05,1.5), "\Omega_m*" )
Obh2_par = Parameter("Obh2", 0.02234, 0.0002, (0.02, 0.025), "\Omega_bh^2")
h_par = Parameter("h", 0.6821, 0.05,(0.4, 1.0), "h")
mnu_par = Parameter("mnu", 0.06, 0.1, (0,1.0), "\Sigma m_{\\nu}")
Nnu_par = Parameter("Nnu", 3.046, 0.5, (3.046,5.046), "N_{\\rm eff}")

Ok_par = Parameter("Ok",0.0, 0.1, (-1.5, 1.5), "\Omega_k")
w_par = Parameter("w",   1.0, 0.1,   (-0.5, 2.0), "w_0")
wa_par = Parameter("wa", 1.0, 0.1,  (-0.5, 2.0), "w_a" )
wb_par = Parameter("wb", 0.7, 0.2,  (-2., 3.0), "w_b")
wc_par = Parameter("wc", 0.7, 0.2,   (-3., 5.0), "w_c")

## this is the prefactor parameter c/rdH0
Pr_par = Parameter("Pr",28.6,4, (5, 70), "c/(H_0r_d)")

## Poly Cosmology Parameters
Om1_par= Parameter("Om1",0.0,0.1,(-3,3), "\Omega_1")
Om2_par= Parameter("Om2",0.0,0.1,(-3,3), "\Omega_2")              

## JordiCDM Cosmology Parameters
q_par=Parameter("q",0.0, 0.2,(0,1), "q")
za_par=Parameter("za",3, 1.0,(2,10), "z_a")
zb_par=Parameter("zb",1, 0.5, (0,2), "z_b")
wd_par=Parameter("wd",-1, 0.5, (-2.0,0.0), "w_d")
Od_par = Parameter("Od",0.0, 0.2, (0.0, 1.0), "O_d")

#Weird model
mu_par = Parameter("mu", 2.13, 0.1, (0.01,5), "\mu")
Amp_par = Parameter("Amp", -0.1839, 0.01, (-2,2), "Amp")
sig_par = Parameter("sig", 0.1, 0.2, (0,3.0), "\sig")

#Tired Light
beta_par = Parameter("beta", 1, 0.1, (-1,1), "\\beta")

#Spline reconstruction
Sp1_par = Parameter("Sp1", 1, 0.01, (-1.2,1.2), "S1")
Sp2_par = Parameter("Sp2", 1, 0.01, (-1.2,1.2), "S2")
Sp3_par = Parameter("Sp3", 1, 0.01, (-1.2,1.2), "S3")
Sp4_par = Parameter("Sp4", 1, 0.01, (-1.2,1.2), "S4")


# Step-wise dark energy density (edit: JG)
# number of redshift boundaries (fixed)
step_nz_par   = Parameter("StepNZ", 3, 1, (0.,10.), "nz")
# redshift of bin boundaries (fixed)
step_z0_par   = Parameter("StepZ0", 0.5, 0.01, (0.,10.), "z0")
step_z1_par   = Parameter("StepZ1", 1.0, 0.01, (0.,10.), "z1")
step_z2_par   = Parameter("StepZ2", 1.6, 0.01, (0.,10.), "z2")
step_z3_par   = None
step_z4_par   = None

# rho_DE/rho_c in redshift bins (nz+1, free), 0.2 is working 
step_rho0_par = Parameter("StepR0", 0.7, 0.5, (-20.,20.), "\\rho_0")
step_rho1_par = Parameter("StepR1", 0.7, 0.5, (-20.,20.), "\\rho_1")
step_rho2_par = Parameter("StepR2", 0.7, 0.5, (-20.,20.), "\\rho_2")
step_rho3_par = Parameter("StepR3", 0.7, 0.5, (-20.,20.), "\\rho_3")
step_rho4_par   = None
step_rho5_par   = None

#Decaying Dark Matter
lambda_par = Parameter("lambda",0.1, 0.2, (0,5.0), "\lambda")
# let the fraction be limited by 0.1 at the bottom.
xfrac_par    = Parameter("xfrac", 1.0, 0.1, (0.1, 1.0), "f_x")

#Early Dark Energy
Ode_par   = Parameter("Ode", 0.05, 0.01, (0,0.9), "\Omega^e_{\\rm de}")

#Slow-Roll DE
dw_par = Parameter("dw",0, 0.1, (-1.0,1.0), "\delta w_0")

#QuinDE
#12, 0.00525, 1,  #4, 0.0625, 1
lam_par = Parameter("lam", 8., 1.5, (0, 20), "\lambda")
V0_par = Parameter("V0", 4.5, 1., (0, 30), "V_0")
A_par  = Parameter("A", 0.00625, 0.005, (0., 0.1), "A")
lB_par 	= Parameter("lB", 1., 0.2, (0, 6), "lB")

#wDark Matter
wDM_par = Parameter("wDM",0.0, 0.1, (-1.0,1.0), "w_{DM}")

#Stiff Cosmology
OSti_par = Parameter("OSti", 0.0, 0.01, (0, 0.1), "O_{Stif}")

#Cosine Parameterisation
nk_par = Parameter("nk", 0.0, 0.1, (-1, 1), "n")
kk_par = Parameter("kk", 1.5, 0.1, (0.01, 3.14), "k")

#Quintom Cosology
mquin_par  = Parameter("mquin", 0.5, 0.1, (0, 3), "m_{\phi}")
mphan_par  = Parameter("mphan", 0.5, 0.1, (0, 3), "m_{\psi}")
iniphi_par = Parameter("iniphi", 0.1, 0.1, (0, 3), "phi_0")

#Anisotropic
bd_par     = Parameter("bd", 1.0, 0.5, (0, 5), "\omega")
Osigma_par  = Parameter("Osigma", -5., 1.0, (-15, 0), "\Omega_{\sigma}")

#Fourier
a0_par     = Parameter("a0", 1, 0.1, (-1, 2), "a0")
a1_par     = Parameter("a1", 0, 0.1, (-1, 2), "a1")
a2_par     = Parameter("a2", 0, 0.1, (-1, 2), "a2")
a3_par     = Parameter("a3", 0, 0.1, (-1, 2), "a3")
a4_par     = Parameter("a4", 0, 0.1, (-1, 2), "a4")

#Horndeski
c1_par     = Parameter("c1", 0.0, 0.15, (-2.0, 1.5), "c1")
c2_par     = Parameter("c2", 0.0, 0.1, (-1.5, 1.5), "c2")
c3_par     = Parameter("c3", 0.0, 0.1, (-1.0, 1.0), "c3")
c4_par     = Parameter("c4", 0.0, 0.1, (-0.5, 0.5), "c4")

f1_par     = Parameter("f1", 1., 0.1, (0.5, 1.5), "f1")
f2_par     = Parameter("f2", 1., 0.1, (0.5, 1.5), "f2")
f3_par     = Parameter("f3", 1., 0.1, (0.5, 1.5), "f3")
f4_par     = Parameter("f4", 1., 0.1, (0.5, 1.5), "f4")





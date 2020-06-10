#!/usr/bin/env python
##
## This file checks against cosmomc results
## AV results are
"""
 -log(Like) =    11.4427384151240
  chi-sq    =    22.8854768302480

   20  0.3080426E+01   logA                  {\rm{ln}}(10^{10} A_s)

    1  0.2230000E-01   omegabh2              \Omega_b h^2
    2  0.1100000E+00   omegach2              \Omega_c h^2
    3  0.1036200E+01   theta                 100\theta_{MC}
    4  0.9000000E-01   tau                   \tau
    5  0.0000000E+00   omegak                \Omega_K
    6  0.6000000E-01   mnu                   \Sigma m_\nu
    7  0.0000000E+00   meffsterile           m_{\nu,{\rm{sterile}}}^{\rm{eff}}
    8 -0.1000000E+01   w                     w
    9  0.0000000E+00   wa                    w_a
   10  0.3046000E+01   nnu                   N_{eff}
   11  0.2400000E+00   yhe                   Y_{He}
   12  0.0000000E+00   alpha1                \alpha_{-1}
   13  0.5000000E+00   deltazrei             {\Delta}z_{re}
   14  0.1000000E+01   Alens                 A_{L}
   15  0.0000000E+00   fdm                   \epsilon_0 f_d
   16  0.1000000E+01   xrs                   Extra
   17  0.9600000E+00   ns                    n_s
   18  0.0000000E+00   nt                    n_t
   19  0.0000000E+00   nrun                  n_{run}
   21  0.0000000E+00   r                     r_{0.05}
   22  0.1000000E+01   Aphiphi               A^{\phi\phi}_L

   23  0.7239477E+00   omegal                \Omega_\Lambda
   24  0.2760523E+00   omegam                \Omega_m
   25  0.0000000E+00   sigma8                \sigma_8
   26  0.1080185E+02   zrei                  z_{re}
   27  0.0000000E+00   r10                   r_{10}
   28  0.6939697E+02   H0                    H_0
   29  0.0000000E+00   r02                   r
   30  0.2176768E+01   A                     10^9 A_s
   31  0.1329451E+00   omegamh2              \Omega_m h^2
   32  0.9225991E-01   omegamh3              \Omega_m h^3
   33  0.2478098E+00   yheused               Y_P
   34  0.1818189E+01   clamp                 10^9 A_s e^{-2\tau}
   35  0.6451439E-03   omeganuh2             \Omega_\nu h^2
   36  0.1389605E+02   age                   {\rm{Age}}/{\rm{Gyr}}
   37  0.1089220E+04   zstar                 z_*
   38  0.1471391E+03   rstar                 r_*
   39  0.1036347E+01   thetastar             \theta_*
   40  0.1059132E+04   zdrag                 z_{\rm{drag}}
   41  0.1498878E+03   rdrag                 r_{\rm{drag}}
   42  0.1378166E+00   kd                    k_D
   43  0.1605557E+00   thetad                \theta_D
   44  0.3161492E+04   zeq                   z_{\rm{eq}}
   45  0.8551706E+00   thetaeq               \theta_{\rm{eq}}
   46  0.1155658E+02   DArs24                D_A(2.4)/(r_s*ext)
   47  0.1180904E+00   H024rs                H(2.4)*r_s*ext
   48  0.3469783E-01   H0rs                  H(0)*r_s*ext
   49  0.2884128E+02   cH0rs                 c/(H_0*r_s*ext)
   50  0.3189206E-01   rsDv24                r_s*ext/D_V(2.4)
   51  0.7748661E-03   H024                  H(2.4)
   52  0.1737575E+04   DA24                  D_A(2.4)

 -log(Like)     chi-sq   data
      3.322      6.644   BAO_Andreu: Lya_Cross_Andreu
      3.482      6.964   BAO_Busca: Lya_Auto_Busca
      4.370      8.740   BAO: DR11CMASS
      0.269      0.537   BAO: DR11LOWZ


> the updated table for the likelihoods is
> 
> -log(Like)     chi-sq   data
>      1.661      3.322   BAO_Andreu: Lya_Cross_Andreu
>      1.741      3.482   BAO_Busca: Lya_Auto_Busca   
>      4.370      8.740   BAO: DR11CMASS
>      0.269      0.537   BAO: DR11LOWZ 
>      

"""
##  
##

from RunBase import *

likelist=[DR11LOWZ(),
    DR11CMASS(),
    DR11LyaAuto(),
    DR11LyaCross()]
likename=['LOWZ','CMASS','LyaAuto','LyaCross']
cosmomc_value=[ 0.269,  4.370,  1.741, 1.661 ]


obh2=0.2230000E-01
om=0.2760523E+00
h=0.6939697
T=LCDMCosmology(obh2,om,h)

print
print
print 'name my_Nlike cosmomc_Nlike, difference'
print '-------------------------------------'
for like,name, cmv in zip(likelist,likename,cosmomc_value):
    like.setTheory(T)
    llike=like.loglike()*-1
    print name, llike, cmv,  llike-cmv



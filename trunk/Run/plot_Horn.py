
from RunBase import *
import numpy as np
import matplotlib.pyplot as plt

T= HornFcoCosmology()

#z = np.linspace(0, 3, 20)
#a = 1./(1+z)
#lna = np.log(a)

eta =  np.linspace(np.log( 1.0/(1.0) ), np.log( 1.0/(4.) ), 20)
a  = np.exp(eta)
z = 1./a -1

print eta
T.ini_f_c()
T.ini_lambda()
h = np.array([T.RHSquared_a(i) for i in a])
T.itype = 'linear'
T.ini_f_c()
T.ini_lambda()
h2 = np.array([T.RHSquared_a(i) for i in a]) #T.hubble(eta[::-1])

dir ='/Users/josevazquezgonzalez/Desktop/Desktop_Jose/work/Papers/Scalar_Fields/Scalar_Fields/'
dataHz = np.loadtxt(dir+'Hz_all.dat')
redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
plt.errorbar(redshifts, obs, errors, xerr=None,
              color='purple', marker='o', ls='None',
              elinewidth =2, capsize=5, capthick = 1, label='$Datos$')


#c = [T.rhow(i) for i in z]
#int_co = [T.Integral_L(i)  for i in z]

#print T.inter(1.0)

#plt.plot(z, int_co)
plt.plot(z, h)
plt.plot(z, h2, label='linear')
#plt.xlim(0, 3.5)
#plt.ylim(0.6, 0.75)
plt.legend(loc = 'best')
plt.show()
#plt.plot(z, h)

#w_= w_par
#wa_=wa_par
#wb_=wb_par
#wc_=wc_par

#w_.setValue(1.6)
#wa_.setValue(0.07)
#wb_.setValue(1.8)
#wc_.setValue(0.0)
#T.updateParams([w_, wa_, wb_, wc_])

#h2 = np.array([np.sqrt(T.RHSquared_a(i))*70 for i in a])
#plt.plot(z, h2, 'r--')
#plt.xlim(0, 3.0)
#plt.ylim(0.6, 300)
#plt.show()

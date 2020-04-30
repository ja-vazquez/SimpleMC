import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt




def init_conditions(a, alpha):
    phi0 = (2 * alpha * (alpha + 2) / 3) ** 0.5 * a ** (3 / (2 + alpha))
    phi1 = 2 * (2 * alpha / ((alpha + 2) * 3)) ** 0.5 * a ** (-3 * alpha / (2 * (2 + alpha)) + 1)
    print phi0, phi1
    return phi0, phi1


def model(Phi, a, H0, Om, Or, alpha, V1cap):
    phi, phiprime = Phi
    t1 = 1 - (a * phiprime) ** 2 / 6
    t2 = V1cap * alpha * phi ** (-(alpha + 1))
    t3 = V1cap * phi ** (-alpha)
    num = t1 * (
    (-2.0 * Or * a ** (-4) - 3 * Om * a ** -3 + t2 / 3 * a * phiprime - 2 * t3) * a * phiprime + 2 * t2 * t3)
    den = 2.0 * (Or * a ** -4 + Om * a ** -3 + t3 / 3) * a ** 2.0
    phidoubleprime = num / den - phiprime / a
    return np.array((phiprime, phidoubleprime))


# def model(Phi, a, H0, Om, Or, alpha, V1cap):
#	phi, phiprime = Phi
#	# to convert to wrt ln(a) to a
#	phiprime *= a
#	t1 = 1 - phiprime**2/6
#	t2 = V1cap*alpha*phi**(-(alpha + 1))
#	t3 = V1cap*phi**(-alpha)
#	num = t1*(t2*t1 - (t3 - Or*a**-4 - t2*phiprime/3)*phiprime)
#	denom = (Om*a**(-3) + Or*a**-4 + t3/3)*(1 + phiprime**2/6)
#	phidoubleprime = (num/denom - phiprime)/a**2
#	return np.array((phiprime, phidoubleprime))

def solvephicdm(H0, Om, Or, alpha, V1cap):
    a_in = 5e-5
    phi0, phiprime0 = init_conditions(a_in, alpha)
    a = np.linspace(a_in, 1, 20000)
    V1capprev = 0
    V1capnew = V1cap
    while abs(V1capnew - V1capprev) > 1e-4:
        V1capprev = V1capnew
        Phi_model = odeint(model, [phi0, phiprime0], a, args=(H0, Om, Or, alpha, V1capnew))
        phi1 = Phi_model[-1, 0]
        phiprime1 = Phi_model[-1, 1]
        V1capnew = 3 * phi1 ** alpha * (1 - Om - Or - phiprime1 ** 2 / 6)
        print V1capnew
    return Phi_model, V1capnew


def getphiandphiprimes(zdata, Phi_model):
    if type(zdata) is float or type(zdata) is np.float64:
        alen = 1
        adata = round(1 / (1 + zdata), 4)
        adata_pos = adata * 20000 - 1
        return Phi_model[int(adata_pos), 0], Phi_model[int(adata_pos), 1]
    else:
        alen = len(zdata)
        adata = 1 / (1 + zdata)
        for i in range(alen):
            adata[i] = round(adata[i], 4)
        adata_pos = adata * 20000 - 1
        phi = np.zeros(alen)
        phiprime = np.zeros(alen)
        for i in range(alen):
            phi[i] = Phi_model[int(adata_pos[i]), 0]
            phiprime[i] = Phi_model[int(adata_pos[i]), 1]
        return phi, phiprime



H0 =70
Om = 0.3
Or =0.0001
alpha = 0.05
V1cap = 1

z = np.linspace(0, 2, 100)
sol = solvephicdm(H0, Om, Or, alpha, V1cap)[0]
phi, phiprime = getphiandphiprimes(z, sol)

a = 1./(1+z)
t1 = 1 - (a * phiprime) ** 2 / 6
t3 = V1cap * phi ** (-alpha)
H2 =H0*np.sqrt((Or * a ** -4 + Om * a ** -3 + t3 / 3)/t1)


plt.plot(z, H2)
#plt.plot(z, phi)
#plt.plot(z, phiprime)
dataHz = np.loadtxt('../data/Hz_all.dat')
redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
plt.errorbar(redshifts, obs, errors, xerr=None,
			 color='purple', marker='o', ls='None',
			 elinewidth =2, capsize=5, capthick = 1, label='$Datos$')
plt.xlabel(r'$z$')
plt.ylabel(r'$H(z) [km/s Mpc^{-1}]$')
plt.title('$\\alpha$ = %0.1f'%alpha)
#plt.savefig('phiCDM.pdf')
plt.show()

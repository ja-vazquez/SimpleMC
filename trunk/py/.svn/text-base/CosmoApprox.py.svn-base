## Various approximations as delivered by M. White
## and rd calculators in the second section
##

import math as N
Tcmb= 2.7255

def z_lastscattering(self,wm,wb):
    """
    z_lastscattering(self,wm,wb):
    Returns z_LS from Hu & White, DampingTail paper.
    """
    b1 = 0.0783*wb**(-0.238)/(1+39.5*wb**0.763)
    b2 = 0.560/(1+21.1*wb**1.81)
    zls= 1048.*(1+0.00124*wb**(-0.738))*(1+b1*wm**b2)
    return(zls)


def z_drag(self,wm,wb):
    """
    z_drag(self,wm,wb):
    Returns z_drag.
    """
    b1 = 0.313*wm**(-0.419)*(1+0.607*wm**0.674)
    b2 = 0.238*wm**0.223
    zd = 1291.*(wm**0.251/(1+0.659*wm**0.828))*(1+b1*wb**b2)
    return(zd)


def soundhorizon_star(self,wm,wb):
    """
    soundhorizon_star(self,wm,wb):
    A fit to the sound horizon, in Mpc, from Eistenstein & Hu (1998;
    ApJ, 496, 605), Eqs. 2-6, except using z_lastscattering not zdrag.
    """
    zeq= 2.50e4 *wm*(self.Tcmb/2.7)**(-4)
    keq= 7.46e-2*wm*(self.Tcmb/2.7)**(-2)      # In 1/Mpc.
    b1 = 0.313*wm**(-0.419)*(1+0.607*wm**0.674)
    b2 = 0.238*wm**0.223
    zd = self.z_lastscattering(wm,wb)
    Rs = 31.5*wb*(self.Tcmb/2.7)**(-4)*(1e3/zd )
    Req= 31.5*wb*(self.Tcmb/2.7)**(-4)*(1e3/zeq)
    s  = 2./3./keq*(6./Req)**0.5*\
         N.log( (N.sqrt(1+Rs)+N.sqrt(Rs+Req))/(1+Req**0.5) )
    return(s)


def soundhorizon_eh(self,wm,wb):
    """
    soundhorizon_eh(self,wm,wb):
    A fit to the sound horizon, in Mpc, from Eistenstein & Hu (1998;
    ApJ, 496, 605), Eqs. 2-6.
    """
    zeq= 2.50e4 *wm*(self.Tcmb/2.7)**(-4)
    keq= 7.46e-2*wm*(self.Tcmb/2.7)**(-2)      # In 1/Mpc.
    b1 = 0.313*wm**(-0.419)*(1+0.607*wm**0.674)
    b2 = 0.238*wm**0.223
    zd = 1291.*(wm**0.251/(1+0.659*wm**0.828))*(1+b1*wb**b2)
    Rd = 31.5*wb*(self.Tcmb/2.7)**(-4)*(1e3/zd )
    Req= 31.5*wb*(self.Tcmb/2.7)**(-4)*(1e3/zeq)
    s  = 2./3./keq*(6./Req)**0.5*\
         N.log( (N.sqrt(1+Rd)+N.sqrt(Rd+Req))/(1+Req**0.5) )
    return(s)


## RD calculators
## These approximation need to have the same interface --
## see how they are used in LCDMCosmology.py

def rd_anderson_approx(obh2,ocbh2,onuh2,Nnu):
    if (abs(Nnu-3)>0.1):
        print "ERROR, cannot use anderson approx with Nnu"
        print "Nnu=",Nnu
    return 55.234 / (ocbh2**0.2538 * obh2**0.1278 * (1+onuh2)**0.3794)

def rd_cuesta_approx(obh2,ocbh2,onuh2,Nnu):
    if (abs(Nnu-3)>0.1):
        print "ERROR, Tony Cuesta says: 'not in this ceral box.'"
        print "Nnu=",Nnu
        stop()
    return 55.154/ (ocbh2**0.25351 * (obh2)**0.12807 * N.exp(
            (onuh2+0.0006)**2.0/0.1176**2 ) )

def rd_cuesta_Nnu_approx(obh2,ocbh2,onuh2,Nnu):
    return 1715.43/(ocbh2)**0.2436/obh2**0.128876/(Nnu-3.046+30.6)/N.exp(49.7*(onuh2+0.002)**2)


def rd_EH_approx(obh2,ocbh2,onuh2,Nnu):
    if (abs(Nnu-3)>0.1):
        print "ERROR, cannot use EH approx with Nnu."
        print "Nnu=",Nnu
        stop()
    return self.CA.soundhorizon_eh(ocbh2, obh2,Nnu)

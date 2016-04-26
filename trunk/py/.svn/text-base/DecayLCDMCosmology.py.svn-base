## This is a CDM cosmology with a decaying
## dark matter component.
##

from pylab import *
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d

from LCDMCosmology import *

class DecayLCDMCosmology(LCDMCosmology):
    ## note that if we don't varyOr, it will be set so that
    ## density at early a is zero.
    def __init__(self, varylam=True, varyxfrac=True, xfrac=xfrac_par.value):
        
        self.varylam=varylam
        self.varyxfrac= varyxfrac
        
        self.lam = lambda_par.value
        self.xfrac = xfrac

        LCDMCosmology.__init__(self)

        self.logar= linspace(0.0,-7.1,100)
        self.ilogar= self.logar[::-1]

        ## force caching
        self.updateParams([])


    ## my free parameters. We add lam, xfrac
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varylam): l.append(lambda_par)
        if (self.varyxfrac): l.append(xfrac_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="lambda":
                self.lam=p.value
            if p.name=="xfrac":
                self.xfrac=p.value
        self.SolveEq()
        ## and updated with relevant rd
        self.setrd(self.rd_func_(self.Obh2, self.Ocbh2_early,self.Omnuh2, self.Nnu()))
        assert(abs(self.RHSquared_a(1.0)-1)<1e-4)
        return True


    def H2_rxrr_a(self,a,rx,rr):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return self.Ocb_std/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Or)+rx/a**3+rr/a**4

    def RHS(self,y,lna):
        ## we are solving rx, so that rhox=rx/a**3 and rhor=rr/a**4
        ##
        a=exp(lna)
        H2=self.H2_rxrr_a(a,y[0],y[1])
        H=sqrt(abs(H2))
        factor = self.lam*y[0]/H
        return array([ - factor, + factor * a])


    def SolveEq(self):
        self.Odm=self.Ocb-self.Obh2/(self.h**2)
        self.Ob=self.Ocb-self.Odm
        self.Ocb_std=self.Ob+self.Odm*(1-self.xfrac)
        self.Odm_dec=self.Odm*self.xfrac
        if (self.lam==0):
            self.Or=0
            yinit= array([self.Odm_dec, self.Or])
            sol =odeint(self.RHS,yinit,self.logar)
        else:
            lowr=0
            highr=0.3
            eps=1e3
            c=0
            while (abs(eps)>1e-7):
                mid=(lowr+highr)/2
                self.Or=mid
                yinit= array([self.Odm_dec, mid])
                sol =odeint(self.RHS,yinit,self.logar)
                self.sol=sol
                eps=sol[-1,1]
                #print "tried: ", mid,eps
                c+=1
                if (eps>0):
                    highr=mid
                else:
                    lowr=mid
                if (c>100):
                    print "Still unconverged after 100 tries."
                    print self.Om, self.Or, self.lam
                    stop()

        ## stupid interp1d doesn't take inverses
        self.sol=sol
        self.rx=interp1d(self.ilogar,sol[::-1,0])
        self.rr=interp1d(self.ilogar,sol[::-1,1])
        ## take early time solution
        self.Ocbh2_early=(self.Ocb_std+sol[-1,0])*self.h**2        
	
    
    def RHSquared_a(self,a):
        lna=log(a)
        return self.H2_rxrr_a(a,self.rx(lna),self.rr(lna))

    def WangWangVec(self):
        print "no WW with Decay"
        stop()
        return None

    ## this returns the "SimpleCMB" variables in a vec
    def CMBSimpleVec(self):
        zstar=1090
        Dastar=self.Da_z(zstar)*self.c_/(self.h*100)
        return array([self.Obh2, self.Ocbh2_early, Dastar/self.rd])
            

    def Om_z(self,a):
	lna=log(a)
	return self.Ocb_std/a**3 + self.rx(lna)/a**3


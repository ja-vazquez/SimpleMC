## Base Cosmology class doesn't know about your
## parameterization of the equation of state or densities or anything.
## However, it does know about Hubble's constant at z=0 OR the prefactor 
## c/(H0*rd) which should be fit for in the case of "rd agnostic" fits.
## That is why you should let it declare those parameterd based on its settings
##
## However, to get the angular diameter distance you need to pass it
## its Curvature parameter (Omega_k basically), so you need to update
## it 


from scipy import *
from RadiationAndNeutrinos import *
from Parameter import *
from ParamDefs import *
import scipy.integrate as intg
from scipy.interpolate import griddata
from scipy.integrate import odeint
import CosmoApprox as CA

class BaseCosmology:
    # speed of light
    c_=299792.45

    def __init__(self,h=h_par.value):
        self.Curv=0
        self.rd=149.50
        self.h=h
        self.prefact=Pr_par.value
        self.varyPrefactor=False
        BaseCosmology.updateParams(self,[])



    def setCurvature(self,R):
        self.Curv=R
    
    def setrd(self,rd):
        self.rd=rd

    def setVaryPrefactor(self,T=True):
        self.varyPrefactor=T


    def setPrefactor(self,p):
        self.prefact=p

    def prefactor(self):
        if self.varyPrefactor:
            return self.prefact
        else:
            return self.c_/(self.rd*self.h*100)

    def freeParameters(self):
        if (self.varyPrefactor):
            Pr_par.setValue(self.prefact)
            l=[Pr_par]
        else:
            h_par.setValue(self.h)
            l=[h_par]
        return l

    def printFreeParameters(self):
        print "Free parameters:"
        self.printParameters(self.freeParameters())

    def printParameters(self,params):
        for p in params:
            print p.name,'=',p.value,'+/-',p.error

    def updateParams(self,pars):
        for p in pars:
            if p.name=="h":
                self.h=p.value
            elif p.name=="Pr":
                self.setPrefactor(p.value)
                ## h shouldn't matter here
                ## we do not want it to enter secondarily through
                ## say neutrinos, so let's keep it 
                ## sane
                #self.h=p.value*self.rd*100/self.c_
        return True
        
    def prior_loglike(self):
        return 0 


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        print "You should not instatiate BaseCosmology"
        error("BAD")
        
    def Hinv_z(self,z):
        return 1/sqrt(self.RHSquared_a(1.0/(1+z)))

    #@autojit
    def DistIntegrand_a(self,a):
        return 1/sqrt(self.RHSquared_a(a))/a**2

    #@autojit
    def Da_z(self,z):
        #r=intg.quad(self.Hinv_z,0,z)
        # this version seems to be faster
        r=intg.quad(self.DistIntegrand_a,1./(1+z),1)

        r=r[0] ## assume precision is ok
        if self.Curv==0:
            return r
        elif (self.Curv>0):
            q=sqrt(self.Curv)
            ## someone check this eq
            ## Pure ADD has a 1+z fact, but have
            ## comoving one
            return sinh(r*q)/(q)
        else:
            q=sqrt(-self.Curv)
            return sin(r*q)/(q)

    # D_a / rd
    def DaOverrd(self,z):
        return self.prefactor()*self.Da_z(z)
    # H^{-1} / rd
    def HIOverrd(self,z):
        return self.prefactor()*self.Hinv_z(z)
    # Dv / rd
    def DVOverrd(self,z):
        return self.prefactor()*(self.Da_z(z)**(2./3.)*(z*self.Hinv_z(z))**(1./3.))

    ## distance modulus
    def distance_modulus(self,z):
        assert(not self.varyPrefactor)
        ### note that our Da_z is comoving, so we're only
        ### multilpyting with a single (1+z) factor
        return 5*log10(self.Da_z(z)*(1+z))

    ## returns the growth factor as a function of redshift

    def GrowthIntegrand_a(self,a):
        return 1/(self.RHSquared_a(a)*a*a)**(1.5)

    def growth(self,z):
        ## Equation 7.80 from Doddie
        af=1/(1.+z)
        r=intg.quad(self.GrowthIntegrand_a,1e-7,af)
        gr=sqrt(self.RHSquared_a(af))*r[0] ## assume precision is ok
        ## If we have Omega_m, let's normalize that way
        if hasattr(self,"Om"):
            gr*=5/2.*self.Om
        return gr



    

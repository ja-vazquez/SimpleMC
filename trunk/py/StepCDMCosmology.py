#
# This is a cosmology with steps in rho_DE(z).
#

import numpy
from numpy import linspace
#from scipy.interpolate import InterpolatedUnivariateSpline
from LCDMCosmology import *
import math as N

class StepCDMCosmology(LCDMCosmology):
    def __init__(self):
      
        self.NZ=int(step_nz_par.value)
        self.Z=numpy.zeros((self.NZ)) # redshifts of bin boundaries 
        self.R=numpy.zeros((self.NZ+1)) # rho_de/rho_c in bins
        
        # see ParamDefs.py for the values of parameters
        if(self.NZ>0) :
            self.Z[0] = step_z0_par.value
            self.R[0] = step_rho0_par.value
            self.R[1] = step_rho1_par.value
        if(self.NZ>1) :
            self.Z[1] = step_z1_par.value
            self.R[2] = step_rho2_par.value
        if(self.NZ>2) :
            self.Z[2] = step_z2_par.value
            self.R[3] = step_rho3_par.value
        if(self.NZ>3) :
            self.Z[3] = step_z3_par.value
            self.R[4] = step_rho4_par.value
        if(self.NZ>4) :
            self.Z[4] = step_z4_par.value
            self.R[5] = step_rho5_par.value
        
        LCDMCosmology.__init__(self)
        
        # compute omegak, updated in updateParams
        self.Ok = 0
        self.Ok = 1.-self.RHSquared_a(1.)
        self.setCurvature(self.Ok)
        
        

        

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if(self.NZ>0) :
            l.append(step_rho0_par)
            l.append(step_rho1_par)
        if(self.NZ>1) :
            l.append(step_rho2_par)
        if(self.NZ>2) :
            l.append(step_rho3_par)
        if(self.NZ>3) :
            l.append(step_rho4_par)
        if(self.NZ>4) :
            l.append(step_rho5_par)
        return l

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            for i in range(self.NZ+1) :
                if p.name=="StepR%d"%i:
                    self.R[i]=p.value
        
        # compute omegak, used in RHSquared_a
        self.Ok = 0
        self.Ok = 1-self.RHSquared_a(1.)
        self.setCurvature(self.Ok)
        return True

    def Rho_de(self,a):
        z=1.0/a-1.0
        for i in range(self.NZ) :
            if (z < self.Z[i]) :
                return self.R[i]
        return self.R[self.NZ]
    
    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        val=self.Ocb/a**3+self.Omrad/a**4+self.Ok/a**2+NuContrib+self.Rho_de(a)
        
        eps=1e-12 # to avoid crazy behavior
        if val<eps :
            val=eps
        return val





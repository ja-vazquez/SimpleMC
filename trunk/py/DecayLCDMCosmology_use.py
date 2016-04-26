## This is LCDM cosmology with optional
## curvature which you can set up with 
## setVaryOk()


from pylab import *
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot  as pyplot
from scipy.interpolate import interp1d

from LCDMCosmology import *

class DecayLCDMCosmology(LCDMCosmology):
    def __init__(self, varyLambda=True, varyOmr=True):
        
        self.varyLambda=varyLambda
        self.varyOmr=varyOmr      
 
        self.Lambda= Lambda_par.value
        self.Omr   =  Omr_par.value

        LCDMCosmology.__init__(self, mnu=0)

        self.solve = self.SolveEq()        

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyLambda): l.append(Lambda_par)
        if (self.varyOmr): l.append(Omr_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="Lambda":
                self.Lambda=p.value
            if p.name=="Omr":
                self.Omr=p.value
        self.solve = self.SolveEq()
        return True


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
#    def RHSquared_a(self,a):         
#        NuContrib=self.NuDensity.rho(a)/self.h**2
#        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*a**(-3*(1.0+self.w)))

#    def Hub2(self,y,lna):
#       return ((1.0-self.Om-self.Omrad)+y[:,0]+y[:,1])

                        #Integrating in log a
    def RHS(self,y,lna):
        factor = self.Lambda*y[0]/(1*self.h*((1.0-self.Om-self.Omr)+ y[0]+ y[1])**(0.5))
        return array([ -3*y[0]/1 - factor, -4*y[1]/1 + factor])


    def SolveEq(self):
       time= linspace(0.0,-7.1,100)
       yinit= array([self.Om, self.Omr])
       self.solve =odeint(self.RHS,yinit,time)
#       print 'hii',self.Lambda, self.Omr
       return self.solve

    

    def Spline(self, lna):         
        y = self.solve        
        x= linspace(-7.1,0.0,100)
        yx= y[:,0] #/self.Hub2(y,x)
        yr= y[:,1] #/self.Hub2(y,x)       #Useful for plotting
        y2x,y2r  =yx[::-1], yr[::-1]
        f = interp1d(x, y2x)
        g = interp1d(x, y2r)
        return f(lna), g(lna)
    


    def RHSquared_a(self,a):
        lna = np.log(a)
        rhox, rhor = self.Spline(lna)
        return (rhox +rhor +(1.0-self.Om-self.Omr))


#Useful for plotting and comparison to LCDM

#    def RHSquared_a2(self,a): 
#        rhox = self.Om/a**3
#        rhor = self.Omr/a**4
#        return (rhox +rhor +(1.0-self.Om-self.Omr))

#    def printA(self):
#       time= logspace(0.0,-5.0,100)
#       y1 = self.SolveEq()
#       fig = pyplot.figure()
#       ax = fig.add_subplot(1,1,1)
#       ax.set_xscale('log')
#       ax.plot(time,y1[:,0]/self.Hub2(y1,time),label='$O_x$')
#       ax.plot(time,y1[:,1]/self.Hub2(y1,time),label='$O_r$')
#       show()

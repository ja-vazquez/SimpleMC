from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import Ochi0_par, r_chi_par, sigma_par, zcr_par

class BECCosmology(LCDMCosmology):
    """
    BEC Cosmology model with Bose-Einstein condensate contribution.
    This model extends the LCDM framework to include a dynamic BEC component for dark matter.
    """
    def __init__(self):       
        # Initialize standard cosmological parameters
        self.Ochi0 = Ochi0_par.value
        self.r_chi = r_chi_par.value
        self.sigma = sigma_par.value
        self.zcr = zcr_par.value
        
        # Call the base class constructor
        LCDMCosmology.__init__(self)

    def freeParameters(self):
        # Add custom parameters for the BEC model to the list of free parameters
        l = LCDMCosmology.freeParameters(self)
        l.append(Ochi0_par)
        l.append(r_chi_par)
        l.append(sigma_par)
        l.append(zcr_par)
        return l

    def updateParams(self, pars):
        # Update standard and custom parameters
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "Ochi0":
                self.Ochi0 = p.value
            if p.name == "r_chi":
                self.r_chi = p.value
            elif p.name == "sigma":
                self.sigma = p.value
            elif p.name == "zcr":
                self.zcr = p.value
        return True

    def RHSquared_a(self, a):
        """
        Modify the Hubble parameter squared H(z)^2/H(z=0)^2, which now includes BEC dark matter.
        """
        bec_term = (self.Ochi0 * self.r_chi/a**3) / ((1 + self.sigma) * (1 + self.zcr)**3) * \
                   (1 - ((self.sigma * a**(-3)) / ((1 + self.sigma) * (1 + self.zcr)**3)))**(-1)
        
        # Hubble parameter squared
        return (self.Obh2/self.h**2)/a**3 + self.Omrad/a**4 + bec_term + (1- self.Ochi0 - self.Omrad )

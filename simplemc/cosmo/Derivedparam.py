

##
# This file has derived parameter definitions given the
# parameters used in this code.
##


class AllDerived:
    """
    Given base parameters, return a list of some Derived ones.
    """
    def __init__(self):

        self.Ol   = Derivedparam('Ol',    0, '\Omega_\Lambda*')
        self.H0   = Derivedparam('H0',    0, 'H_0*')
        self.Age  = Derivedparam('Age',   0, 'Age[Gyr]*')

        # Parameters for the Omh2 test, to see deviations from
        # LCDM at three redshifts.
        self.Omh12= Derivedparam('Omh12', 0, 'Omh2(z1;z2)*')
        self.Omh13= Derivedparam('Omh13', 0, 'Omh2(z1;z3)*')
        self.Omh23= Derivedparam('Omh23', 0, 'Omh2(z2;z3)*')

        # Parameter from DGP models.
        self.Orc  = Derivedparam('Orc',   0, '\Omega_{rc}*')
        self.list = [self.Ol, self.H0, self.Age, self.Omh12, self.Omh13, self.Omh23, self.Orc]


    def listDerived(self, like):
        """
        Given the free parameters compute derived ones.

        Parameters
        ----------
        like: object
            object defined in BaseLikelihood, that contains
            free-parameters and theory.

        Returns
        -------
        list: list
            List with values from derived parameters.
        """
        self.like  = like
        self.cpars = like.freeParameters()
        self.Ol.setValue(   self.computeDerived('Ol'))
        self.H0.setValue(   self.computeDerived('H0'))
        self.Age.setValue(  self.computeDerived('Age'))
        self.Omh12.setValue(self.computeDerived('Omh12'))
        self.Omh13.setValue(self.computeDerived('Omh13'))
        self.Omh23.setValue(self.computeDerived('Omh23'))
        self.Orc.setValue(  self.computeDerived('Orc'))
        return self.list


    def computeDerived(self, parname):
        """Initialize and compute Derived parameters in
        terms of the base ones."""
        if parname == 'Ol':
            for par in self.cpars:
                if par.name == 'Om': return 1- par.value
        if parname == 'Orc':
            for par in self.cpars:
                if par.name == 'Om': return (1- par.value)**2/4.
        elif parname == 'H0':
            for par in self.cpars:
                if par.name == 'h': return par.value*100
        elif parname == 'Omh12':
            return self.Omh2(0.0, 0.57)
        elif parname == 'Omh13':
            return self.Omh2(0.0, 2.34)
        elif parname == 'Omh23':
            return self.Omh2(0.57, 2.34)
        elif parname == 'Age':
            return self.like.theory_.Age()
        else:
            import sys
            sys.exit('Define derived parameter', parname)




    def Omh2(self, z1, z2):
        """
        Computes the Omh2 diagnostic at two redshifts
        to test deviations from the LCDM model.
        For LCDM the value is a constant for any combination
        of redshifts.

        Parameters
        ----------
        z1: float
            Redshift z1 to compute the diagnotic.

        z2: float
            Redshift z2 to compute the diagnotic.

        Returns
        -------
        Obh2: float
            Obh2 dianotic: see [arXiv:1406.2209].
        """
        h0 = self.like.theory_.h
        h1 = h0**2*self.like.theory_.RHSquared_a(1.0/(1+z1))
        h2 = h0**2*self.like.theory_.RHSquared_a(1.0/(1+z2))

        Obh2 = (h1-h2)/((1+z1)**3 - (1+z2)**3)
        return Obh2



class Derivedparam:
    """
    Auxiliary class, based on Parameter class.

    Parameters
    ----------
    name: string
        Parameter name.

    value: float
        Initialize the value, and use a function to update it.

    Ltxname: string, optional
        Provide the Latex name, useful for plotting.
        Default is None, and in this case uses the 'name' string.

    Example
    -------
    Hubble parameter
    self.H0   = Derivedparam('H0',  function, 'H_0*')
    """
    def __init__(self, name, value, Ltxname=None):

        self.name = name
        if Ltxname:
            self.Ltxname = Ltxname
        else:
            self.Ltxname = name
        self.value = value


    def setLatexName(self, Ltx):
        self.Ltxname = Ltx

    def setValue(self, val):
        self.value = val




class AllDerived:
    def __init__(self):
        """
        Given base parameters, add some Derived ones.
        Returns
        -------

        """
        self.Ol   = Derivedparam('Ol',    0, '\Omega_\Lambda*')
        self.H0   = Derivedparam('H0',    0, 'H_0*')
        self.Age  = Derivedparam('Age',   0, 'Age[Gyr]*')
        self.Omh12= Derivedparam('Omh12', 0, 'Omh2(z1;z2)*')
        self.Omh13= Derivedparam('Omh13', 0, 'Omh2(z1;z3)*')
        self.Omh23= Derivedparam('Omh23', 0, 'Omh2(z2;z3)*')
        self.Orc  = Derivedparam('Orc',   0, '\Omega_{rc}*')
        self.list = [self.Ol, self.H0, self.Age, self.Omh12, self.Omh13, self.Omh23, self.Orc]


    def listDerived(self, like):
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
        h0 = self.like.theory_.h
        h1 = h0**2*self.like.theory_.RHSquared_a(1.0/(1+z1))
        h2 = h0**2*self.like.theory_.RHSquared_a(1.0/(1+z2))
        return (h1-h2)/((1+z1)**3 - (1+z2)**3)



class Derivedparam:
    def __init__(self, name, value, Ltxname=None):
        """
        Auxiliary class, based on Parameter class
        Parameters
        ----------
        name: parameter name
        value: standard value
        Ltxname: Latex name

        Returns
        -------

        """
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
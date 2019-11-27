##
# This is a DV likelihood that comes in form of a table.
##

import scipy as sp
from BaseLikelihood import BaseLikelihood
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import quad


class TabulatedBAODVLikelihood(BaseLikelihood):
    def __init__(self, name, filename, fid_theory, z):
        BaseLikelihood.__init__(self, name)
        print()
        print("Loading ", filename)
        data = sp.loadtxt(filename)
        self.chi2i = interp1d(data[:, 0], data[:, 1])
        self.fidDVOverrd = fid_theory.DVOverrd(z)

        # find minimum
        alphamin = minimize(self.chi2i, 1.0, bounds=[
                            [0.9, 1.1]], method='SLSQP').x[0]
        rms  = quad(lambda x: sp.exp(-self.chi2i(alphamin+x)/2)*x*x, -0.1, 0.1)[0]
        rms /= quad(lambda x: sp.exp(-self.chi2i(alphamin+x)/2), -0.1, 0.1)[0]
        DVf  = self.fidDVOverrd*fid_theory.rd
        print(name, "measurement in", fid_theory.rd_approx,
              " : DV=", DVf*alphamin, "+-", sp.sqrt(rms)*DVf)
        print("with rd=", fid_theory.rd, "DV_fid=", DVf, "alphamin=", alphamin)
        self.z = z


    def loglike(self):
        alpha = self.theory_.DVOverrd(self.z)/self.fidDVOverrd
        try:
            chi2 = self.chi2i(alpha)
        except:
            # print "Note: alpha for ",self.name(),"out of lookup-table bounds"
            chi2 = 9
        return -chi2/2.0

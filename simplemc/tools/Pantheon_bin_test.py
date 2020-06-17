
import numpy as np
import matplotlib.pyplot as plt
from simplemc.tools.cosmich import cosmochain

class GenericPantheona():
    def __init__(self):
        """
        This class is useful to compare compress Pantheon dataset (from binning)
        with the full Phanteon , i.e. see models/CompressPantheon.py
        Returns
        -------

        """

        self.zini    = 0 # 13.9

        nbins = 15

        dir_name   = 'chains/'
        root = 'GPantheon_phy_CPantheon_'+'%d'%nbins+'_mcmc'

        C= cosmochain(dir_name + root)
        params= ['zbin%d'%i for i in range(nbins)]
        x =C.GetCovariance(params)
        self.zs= x[0]
        self.cov= x[1]
        self.best_fit = [C.GetLimits(p, returnlims=True)[3] for p in params]

        self.parvals = np.ones(len(self.zs))
        self.zvals = np.logspace(np.log10(0.01),np.log10(2.261), len(self.zs)+1)
        print('zvals', self.zvals )
        da=[x.split() for x in open('data/pantheon_lcparam_full_long_zhel.txt').readlines()[1:]]
        self.zcmb = np.array([float(line[1]) for line in da])
        self.mag = np.array([float(line[4]) for line in da])



    def mu_bar(self,n,z,m1,m2):
        alpha = np.log10(z/self.zvals[n])/np.log10(self.zvals[n+1]/self.zvals[n])
        return (1- alpha)*m1 + alpha*m2



    def genericPModel(self, z):
        if z>=self.zvals[0] and z< self.zvals[len(self.parvals)]:
            for i in range(len(self.parvals)):
                if z>=self.zvals[i] and z< self.zvals[i+1]:
                    if i ==0:
                        y = self.mu_bar(0,  z, self.zini, self.zs[0])
                    else:
                        y = self.mu_bar(i,  z, self.zs[i-1], self.zs[i])
        else:
            y = 0
        return y


    def plots(self):
        """
        Compare full pantheon with generic model (100 points) and from
        best-fit values.
        Returns
        -------

        """
        plt.plot(self.zcmb, self.mag, 'ro')
        z= [i for i in np.linspace(0.01, 2.26, 100)]
        plt.plot(z, [self.genericPModel(i)+13.9 for i in z])

        self.zs = self.best_fit
        plt.plot(z, [self.genericPModel(i)+13.9 for i in z], 'r')
        plt.show()



    def wfiles(self):
        with open('binned_pantheon_%d.txt'%len(self.zs), 'w') as f:
            f.write('# binned Pantheon SNe\n')
            f.write('# z mu\n')
            for i, j in enumerate(self.zs):
                f.write('%f\t %f \n'%(self.zvals[i+1], j+ 13.9))
        #with open('binned_cov_pantheon.tx', 'w') as g:
        #    g.write('%d\t %d\n'%(len(self.zs), len(self.zs)))
        np.savetxt('binned_cov_pantheon_%d.txt'%len(self.zs), self.cov, header='%d\t %d'%(len(self.zs), len(self.zs)))


if __name__=='__main__':
    A= GenericPantheona()
    A.plots()
    #A.wfiles()
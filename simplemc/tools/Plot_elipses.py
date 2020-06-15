import scipy as sp
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np

def plot_elipses(best, cov, par1, par2, ax=None, **kwargs):
            #fig = plt.figure(figsize=(6,6))
            #ax = fig.add_subplot(111) #, aspect='equal')

            #sigmas
            sigmas = [2.3, 6.18] #, 11.83]

            #par1, par2 index of the cov matrix
            covslice = cov[[par1,par2], :][:, [par1,par2]]
            mn  = best
            ax  = ax

            def eigsorted(cov):
                vals, vecs = sp.linalg.eigh(cov)
                order = vals.argsort()[::-1]
                return vals[order], vecs[:,order]

            vals, vecs = eigsorted(covslice)
            vec = vecs.T.copy()

            #plot vectors
            vec[0]*=sp.sqrt(6.18*sp.real(vals[0]))
            vec[1]*=sp.sqrt(6.18*sp.real(vals[1]))

            plt.plot(mn[par1], mn[par2],'bo') #, label=self.model)
            plt.plot([mn[par1]-vec[0][0], mn[par1]+vec[0][0]],
                [mn[par2]-vec[0][1],mn[par2]+vec[0][1]],'r-')
            plt.plot([mn[par1]-vec[1][0],mn[par1]+vec[1][0]],
                [mn[par2]-vec[1][1],mn[par2]+vec[1][1]],'r-')

            theta = sp.degrees(np.arctan2(*vecs[:,0][::-1]))

            for i, sigs in enumerate(sigmas):
                w, h = 2*sp.sqrt(vals)*sp.sqrt(sigs)
                ell = Ellipse(xy=(mn[par1], mn[par2]),  width = w, height = h,\
                              angle=theta, color='green',  lw=4)
                ell.set_facecolor('none')
                ax.add_artist(ell)

            ax.legend([ell], ['Fisher'])
            #plt.legend(loc='best')
            #plt.title('Fisher', fontsize=10)
            #plt.show()

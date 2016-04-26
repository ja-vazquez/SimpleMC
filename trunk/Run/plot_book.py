#!/usr/bin/env python

from scipy import *
from cosmich import *
from plot_aux import *
import pylab, sys




params1 = {'backend': 'pdf',
               'axes.labelsize': 32,
               'text.fontsize': 20,
               'xtick.labelsize': 24,
               'ytick.labelsize': 24,
             #  'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)

#-----------
dire = 'chains/'
show=True 


T= plot_aux()



if True:
        datasetl = ['SN','Planck_15','Planck_15+SN','Planck_15+SN+BBAO']
        model    = ['oLCDM']
        extra    = 'phy'
        params   =  ['Om','Ol']
        lw=[4,4,4]
        lcolor=['blue','red','green','yellow']
        T.Plotting_2d(dire, model, extra, datasetl, params, loc='upper left', Ol='True',legcolor=True,Nb=40)
        pylab.ylabel('$\\Omega_{\Lambda}$')
        #pylab.plot([0,1],[1,0],'k:')
        #pylab.xlim(xmin=0.06)
	#pylab.ylim(ymax=1)
        pylab.tight_layout()
        # pylab.savefig('plot_book.pdf')
        if show: pylab.show()

if False:
	m = np.loadtxt("book.txt")
	x = [130.767528128, 2562.40366829, 12570.0411206]
	s = (0.042064972, 0.851764228, -0.167495872, 1.434176672)
	pylab.contourf(m,x,extent=s, origin='lower',aspect='auto', cmap=pylab.get_cmap('Greens'))
        pylab.contour(m,x,extent=s, origin='lower',aspect='auto', colors='green', lwidth=[4,4,4])

	pylab.axis([0.2,0.4,0.5,0.85])
	pylab.show()







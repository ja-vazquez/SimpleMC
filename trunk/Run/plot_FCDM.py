import numpy as np
#from RunBase import *
#from cosmich import *
from ChainIterator import *


if(True):                  #Compute functions= ratio of Fig 2 at z=0.57/0
          dire2 = '/astro/u/jvazquez/work/SimpleMC/trunk/chains/'
          D=ChainIterator(dire2, 'FCDM','phy', 'Planck_15+BBAO+SN', balance =False)
          grlist=[]
          wlist=[]
	  with open('my_data.txt', 'w') as g:
           for i in range(0,D.N,2000):
            T=D.theory(i)
	    
            z = [j for j in np.arange(-0.9,4,0.2)]
	    wde = [T.w_de(j) for j in np.arange(-0.9,4,0.2)]
	    ll = len(z)

            t = (D.weight(i),) + tuple(z) + tuple(wde)
	    g.write(str('%2.5e\t'*(2*ll+1)%(t))) 
	    g.write('\n')



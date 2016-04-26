#conpute convergence values for the MCMC chains
#!/usr/bin/env python
import pylab
from scipy import *
import scipy.linalg as la

Tmean=[0.3084178E+00,0.2252339E-01 ,0.6708606E+00,-0.9510208E+00,-0.8414369E-03]

mean=0
cova=0
weight=0
weightsq=0
b =[]
i=0


dire='chains_SimpleMC/'
file_name='owCDM_phy_BBAO+Planck_3.txt'
outfile='MCMC_owCDM'


formn ='%g '*(2) + '\n'
fpar=open(outfile+".converge",'w')


with open(dire+file_name) as inf:
    for line in inf:
        parts = [float(x) for x in line.split()]
        if len(parts) > 1:   # if at least 2 parts/columns
	     we=parts[0]
             vec= parts[2:7]   
	     	
	     weight+= we
	     weightsq+=we**2
	     mean+= (we*array(vec))
	     a=array(mean/weight)	
	     #b.append(a)

	     cova+= we*array((a-vec)*array([a-vec]).T)		     
	     if i>0:	
	      factors=(weight/(weight**2-weightsq))	

	      fcova = factors*cova 
	      conver = dot(dot(a-Tmean, la.inv(fcova)),a-Tmean) 

	      if(not i%1000):	
                cvg=formn%tuple([i, conver])
                fpar.write(cvg)
                 	     
	     i+=1

fpar.close()

	

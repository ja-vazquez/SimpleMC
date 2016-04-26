#!/usr/bin/env python

# root h2list.C

import sys,string,os.path
import numpy as N
import pylab as P
from StringIO import StringIO   # StringIO behaves like a file object
from subprocess import call
from matplotlib.pyplot import *
import pylab

#convert=True

def mycontour(listfilename,color='g',label='toto',**keys) :
    
    print "reading",listfilename
    
    if not os.path.isfile(listfilename) :
        print "missing",listfilename
        sys.exit(12)
    
    z,param_1sig_minus,param_best,param_1sig_plus = N.loadtxt(listfilename, delimiter=' ', usecols=(0,1,2,3), unpack=True)
    
    colorthemes = {'g': [(0, 0.33, 0),(0, 0.67, 0),(144 / 255., 238 / 255., 144 / 255.)],
                   'b': [(0, 0, 0.8),(0, 0.4, 1.0),(0, 183 / 255., 1)],
                   'k': ['0.2','0.4','0.6'],
                   'r': [(0.3,0,0),(0.6,0.,0), (0.9, 0, 0)],
                   }
    
    colors=colorthemes[color]
    P.plot(z,param_1sig_minus,color=colors[0])
    P.plot(z,param_1sig_plus,color=colors[0])
    res=P.fill_between(z,param_1sig_minus,param_1sig_plus,color=colors[1],alpha=0.5,**keys)
    return res





nodes=[0.,0.11,0.5,2.2]

if len(sys.argv)<2 :
    print sys.argv[0],"plot.conf"
    sys.exit(12)

filenames={}
labels={}
conffilename=sys.argv[1]
file=open(conffilename)
for line in file.readlines() :
    if line[0]=="#" :
        continue
    vals=string.split(line.strip()," ")
    if len(vals)<3 :
        continue
    
    what=vals[0]
    if not filenames.has_key(what) :
        filenames[what]=[]
    if not labels.has_key(what) :
        labels[what]=[]
    
    filenames[what].append(vals[1])
    labels[what].append(vals[2].replace("_"," "))
    

if filenames.has_key("R") : 
    fig=P.figure("dark_energy_density")
    colors=["g","b","r","k"]

    #ylim=[-0.3,1.8]
    ylim=[-1.2,1.2]
    for n in nodes :
        P.plot([n,n],ylim,ls="--",c="Gray",lw=1)
    
    contours=[]
    c=0
    for filename in filenames["R"] :
        print filename,labels["R"][c]
        contours.append(mycontour(filename,color=colors[c],label=labels["R"][c]))
        c+=1
    proxy = [P.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) for pc in contours]

    P.legend(proxy,labels["R"],loc=4,prop={"size":12})
    
    P.ylim(ylim)
    
    

    P.xlabel(r"$z$")
    P.ylabel(r"$\rho_{DE}/\rho_c$")
    suffix=string.replace(conffilename,".conf","")
    suffix=string.replace(suffix,"plot-","")
    figfilename="dark_energy_density_%s.pdf"%suffix
    fig.savefig(figfilename)
    print "wrote",figfilename

if filenames.has_key("H")  : 
    fig=P.figure("inverse_distance_ladder")
    colors=["g","b","r","k"]
    
    ylim = [50,90]
    if False :
        for n in nodes :
            P.plot([n,n],ylim,ls="--",c="Gray",lw=1)

    contours=[]
    c=0
    for filename in filenames["H"] :
        contours.append(mycontour(filename,color=colors[c],label=labels["H"][c]))
        c+=1

    marker_size=10.
    # H0
    z=N.array([0.,0.])
    hz=N.array([73.8,74.3])
    hze=N.array([2.4,2.3])  
    P.errorbar(z,hz,hze,color='k',fmt='o',ms=marker_size,lw=2,label=r"local $H_0$")
    
    # true H(z) gal
    z=N.array([0.57])
    hz=N.array([96.8/1.57]) #?
    hze=N.array([3.4/1.57])  
    P.errorbar(z,hz,hze,color='r',fmt='o',ms=marker_size,lw=2,label=r"$H(z=0.57)$ gal.")
    
    # eff. H(z) gal
    z=N.array([0.204,0.255])
    hz=N.array([60.43,60.857])
    hze=N.array([1.25,0.882])  
    P.errorbar(z,hz,hze,color='r',fmt='^',ms=marker_size,lw=2,label=r"$D_V(z=0.32)$ $D_A(z=0.57)$ gal.")
    
    # eff. H(z) gal 6dF
    z=N.array([0.0682])
    hz=N.array([63.75])
    hze=N.array([2.85])  
    P.errorbar(z,hz,hze,color='g',fmt='^',ms=marker_size,lw=2,label=r"$D_V(z=0.106)$ 6dF")
    
    # eff. H(z) gal WiggleZ
    z=N.array([0.276,0.4247])
    hz=N.array([60.07,64.48])
    hze=N.array([5.46,4.11])  
    #P.errorbar(z,hz,hze,color='b',fmt='^',ms=marker_size,lw=2,label="pseudo. H(z) WiggleZ")
    
    # true H(z) Lyman-alpha
    z=N.array([2.34,2.36])
    hz=N.array([221.6/3.34,226/3.36])
    hze=N.array([6.7/3.34,8/3.36])  
    P.errorbar(z,hz,hze,color='y',fmt='o',ms=marker_size,lw=2,label=r"$H(z=2.34)$ $H(z=2.36)$ Ly-$\alpha$")

    # eff. H(z) Lyman-alpha (see note)
    z=N.array([0.696,0.698])
    hz=N.array([61.57,64.41])
    hze=N.array([3.12,2.10])  
    P.errorbar(z,hz,hze,color='y',fmt='^',ms=marker_size,lw=2,label=r"$D_A(z=2.34)$ $D_A(z=2.36)$ Ly-$\alpha$")
    
    

    proxy = [P.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) for pc in contours]
    l1=P.legend(proxy,labels["H"],loc=4,prop={"size":12})
    l2=P.legend(numpoints=1,loc=9,prop={"size":12})
    gca().add_artist(l1)
    
    P.xlim([-0.1,3])
    P.ylim(ylim)
    
    P.xlabel(r"$z$")
    P.ylabel(r"$H(z)/(1+z)$")
    suffix=string.replace(conffilename,".conf","")
    suffix=string.replace(suffix,"plot-","")
    figfilename="inverse_distance_ladder_%s.pdf"%suffix
    fig.savefig(figfilename)
    print "wrote",figfilename

P.show()

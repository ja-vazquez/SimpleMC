#!/usr/bin/env python

from cosmich import *
PLAdir='/astro/u/anze/work/Planck/PLA'

def get_chain(model, combination):
    return cosmochain("%s/%s/%s/%s_%s"%(PLAdir,model,combination,model,combination))


c1=get_chain('base_omegak','planck_lowl_lowLike')
c2=get_chain('base_nrun_r_omegak','planck_lowl_lowLike')
c3=get_chain('base','planck_lowl_lowLike')

hack=True
if hack:
    okerr=0.02
    oki=c1.parcol['omegak']
    c1.chain[:,0]*=exp(-c1.chain[:,oki]**2/(2*okerr**2))
    oki=c2.parcol['omegak']
    c2.chain[:,0]*=exp(-c2.chain[:,oki]**2/(2*okerr**2))


pylab.figure()

for i,p in enumerate(['omegabh2','omegach2','thetastar*']):
    pylab.subplot (3,2,2*i+1)
    for c,t in zip([c1,c2,c3],['r-','b-','g-']):
        x,y=c.GetHisto(p)
        pylab.locator_params(nbins=6)        
        pylab.plot(x,y,t)
        pylab.locator_params(nbins=6)        
    pylab.subplot(3,2,2*i+2)
    for c,t in zip([c1,c2],[-1,-2]):
        c.Plot2D(p,'omegak', filled=t)

if hack:
    pylab.savefig('pars_hacked.pdf')
else:
    pylab.savefig('pars.pdf')





#!/usr/bin/env python
from RunBase import *
from cosmich import *
PLAdir='/astro/u/anze/work/Planck/PLA'

## model, pick your poison

model='base_omegak'
combination='WMAP'

#model='base_nrun_r_omegak'
#combination='planck_lowl_lowLike'

#model='base_mnu_omegak'

#pick the second poison
# plack=camspec 50<l<2500
# lowl=low planck (l<50 TT)
# lowLike= WMAP low polarization
# Antony, thank you for this amazing clarity




C=cosmochain("%s/%s/%s/%s_%s"%(PLAdir,model,combination,model,combination))
okerr=0.04
oki=C.parcol['omegak']
C.chain[:,0]*=exp(-C.chain[:,oki]**2/(2*okerr**2))


obh2i, och2i, oki, omi,H0i, zstari, rstari, tstari, zdragi, rdragi, kdi, thetadi, zeqi, thetaeqi=[
C.parcol[name] for name in ['omegabh2','omegach2','omegak','omegam*','H0*', 'zstar*', 'rstar*', 'thetastar*','zdrag*', 
'rdrag*', 'kd*', 'thetad*', 'zeq*', 'thetaeq*']]

T=oLCDMCosmology()

mvec=zeros(3)
cmvec=zeros(3)
mmat=zeros((3,3))
cmmat=zeros((3,3))
sw=0

omh2a=0
omh2aa=0


rdrat=[]
N=len(C.chain)

for i,sample in enumerate(C.chain[:]):
    #surely there must be some less shitty way of setting this up.
    obh2, och2, ok, om,H0, zstar, rstar, tstar, zdrag, rdrag, kd, thetad, zeq, thetaeq = [
        sample[ndx] for ndx in [obh2i, och2i, oki, omi,H0i, zstari, rstari, tstari, zdragi, 
        rdragi, kdi, thetadi, zeqi, thetaeqi]]
    #print , om***2,'AA'
    myom=(obh2+och2)/(H0/100.0)**2
    T.updateParams([Parameter("Ok",ok), Parameter("h",H0/100), Parameter("Obh2",obh2),
                    Parameter("Om",myom)])
    rdrat.append(T.rd/rdrag)
    ##our wang wang vec
    vec=T.WangWangVec()
    #recompute from their quantities
    Dastar=100*rstar/tstar
    R=sqrt(om)*H0*Dastar/T.c_
    la=pi*Dastar/rstar
    cvec=array([la,R,obh2])
    # weight
    w=sample[0]
    mvec+=w*vec
    mmat+=w*outer(vec,vec)
    cmvec+=w*cvec
    cmmat+=w*outer(cvec,cvec)
    omh2=om*(H0**2)/10000.
    omh2a+=w*omh2
    omh2aa+=w*omh2**2
    sw+=w

    if (i%100==0):
        print 1.0*i/N*100,"%"

def print_stat(vec,mat):
    vec/=sw
    mat/=sw
    mat-=outer(vec,vec)
    err=sqrt(mat.diagonal())
    print vec
    print err
    print mat
    print mat*outer(1/err,1/err)

print "My quantities:"
print_stat(mvec,mmat)

print "Cosmoc quantities:"
print_stat(cmvec, cmmat)

print "mean Omh2"
omh2a/=sw
omh2aa/=sw
omh2aa-=omh2a**2
print omh2a, sqrt(omh2aa)




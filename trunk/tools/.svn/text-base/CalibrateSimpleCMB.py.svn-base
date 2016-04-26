#!/usr/bin/env python
from RunBase import *
from cosmich import *
PLAdir='/astro/u/anze/work/Planck/PLA'

## model, pick your poison

model='base_Alens'
combination='WMAP'
#combination='planck_lowl_lowLike'


#pick the second poison
# plack=camspec 50<l<2500
# lowl=low planck (l<50 TT)
# lowLike= WMAP low polarization
# Antony, thank you for this amazing clarity




C=cosmochain("%s/%s/%s/%s_%s"%(PLAdir,model,combination,model,combination),'auto')

obh2i, och2i, omi,H0i, zstari, rstari, tstari, zdragi, rdragi, kdi, thetadi, zeqi, thetaeqi=[
C.parcol[name] for name in ['omegabh2','omegach2','omegam*','H0*', 'zstar*', 'rstar*', 'thetastar*','zdrag*', 
'rdrag*', 'kd*', 'thetad*', 'zeq*', 'thetaeq*']]

T=LCDMCosmology()

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
    obh2, och2, H0 = [sample[ndx] for ndx in [obh2i, och2i, H0i]]

    #print , om***2,'AA'
    myom=(obh2+och2)/(H0/100.0)**2
    T.updateParams([Parameter("h",H0/100), Parameter("Obh2",obh2),
                    Parameter("Om",myom)])
    ##our wang wang vec
    vec=T.CMBSimpleVec()
    # weight
    w=sample[0]
    mvec+=w*vec
    mmat+=w*outer(vec,vec)
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




